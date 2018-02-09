import matplotlib
matplotlib.use('TkAgg')

import glob
import argparse
import ephem
import sys
import airballoon
#import time
import numpy as np
from joblib import Parallel, delayed
from photutils import DAOStarFinder, CircularAperture, CircularAnnulus, aperture_photometry, SExtractorBackground
from astropy.io import fits
from astropy.utils.console import color_print, ProgressBar
from astropy.stats import sigma_clipped_stats, SigmaClip
from utils import altaz2xy
from locations import locations

parser = argparse.ArgumentParser(description='Extinction map')
parser.add_argument('Folder', type=str, help='Folder containing .fits files for one night')
parser.add_argument('--location', type=str, default='Gemini', help='Defined location on locations.py')
parser.add_argument('--manual', action='store_true', help='Select manually one star')
parser.add_argument('--star', default='Sirius', type=str, help='Star to to follow (exp)')
parser.add_argument('--follow', action='store_true', help='Show star paths')

args = parser.parse_args()

args.star = args.star.split(',')
loc = locations[args.location]()
folder_date = args.Folder.split('/')[-2]

#Leer archivos
color_print('\nLeyendo archivos...', 'cyan')
print '\tBuscando archivos fits'
all_fits = np.sort(glob.glob('%s*.fits*' % args.Folder))
print '\tEncontrados %d archivos' % len(all_fits)
print '\tLeyendo headers...'
all_head = Parallel(n_jobs=4, verbose=5)(delayed(fits.getheader)(f,1) for f in all_fits)
dates    = np.array([ephem.Date(h['DATE']+ ' ' + h['UT']) for h in all_head])
loc.date = dates[0]

print '\tFiltrando imagenes nocturnas'
loc.horizon = '-12'  #Astronomical Twilight
sunset  = loc.next_setting(ephem.Sun())
sunrise = loc.next_rising(ephem.Sun())
noche   = (dates > sunset) * (dates < sunrise)
if noche.sum==0:
    color_print('\nNo se encontraron imagenes nocturnas', 'cyan')
    sys.exit()

print '\t\tTwilight (Sunset) : ' + str(sunset)
print '\t\tTwilight (Sunrise): ' + str(sunrise)
print '\t\tSe usaran %d imagenes' % noche.sum()

noche[:10]  = False
noche[-10:] = False

all_fits = all_fits[noche]
numeros  = np.arange(len(all_head))[noche]
exptimes = np.array([h['EXPTIME'] for h in all_head])[noche]

print '\t\tPrimera imagen: ' + all_fits[0]
print '\t\tUltima imagen:  ' + all_fits[-1]

print '\tLeyendo datos...'
all_data = np.array(Parallel(n_jobs=4, verbose=5)(delayed(fits.getdata)(f) for f in all_fits))#np.array([fits.getdata(f) for f in ProgressBar(all_fits)])#np.array(Parallel(n_jobs=4, verbose=0)(delayed(fits.getdata)(f) for f in all_fits))

if args.manual:
    import matplotlib.pyplot as plt
    xco = []
    yco = []
    def onclick(event):
        xco.append(event.xdata)
        yco.append(event.ydata)

    flux = np.zeros(noche.sum())
#xco  = np.zeros(noche.sum())
#yco  = np.zeros(noche.sum())

    few_data = all_data[::5]

    fig, ax = plt.subplots()
    for i in range(len(few_data)):
        fig.canvas.mpl_connect('buttom_press_event', onclick) 
        ax.matshow(few_data[i], cmap='gray_r')
        fig.show()

        print xco
        print yco
    sys.exit(1)
    
    

print '\t\tCalculando posiciones de estrellas'
if len(args.star) > 1:
    star = ephem.readdb("%s,f|S|A0,%s,%s,2.00,2000" % tuple(args.star))
else:
    star = ephem.star(args.star[0])
x, y = np.zeros((2, noche.sum()))
flux    = np.zeros(noche.sum())
airmass = np.zeros(noche.sum())
dateflo = np.zeros(noche.sum())
datetim = []

YY, XX = np.ogrid[:all_data[0].shape[0], :all_data[0].shape[1]]
sigma_clip = SigmaClip(sigma=3)

for i in ProgressBar(noche.sum()):
    loc.date = dates[noche][i]
    star.compute(loc)
    el = float(repr(star.alt))#*180/np.pi
    az = float(repr(star.az))#*180/np.pi

    x[i], y[i] = altaz2xy(el, az)

    #Background y unidades
    nowdata = all_data[i]
    nowhead = fits.getheader(all_fits[i], 1)
 
    pedestal = nowhead['PEDESTAL']
    egain    = float(nowhead['EGAIN'])
    exptime  = float(nowhead['EXPTIME'])

    nowdata = nowdata + pedestal
    nowdata[nowdata < 0] = 0
    nowdata = nowdata / exptime

    #Background
    bdist = np.sqrt((XX - x[i])**2 + (YY - y[i])**2)
    bmask = (bdist < 10.)# * (bdist > 4.)

    for j in range(30):
        btest = (bdist < 6+j)
        mtest, vtest, _ = sigma_clipped_stats(nowdata[btest], sigma=3.0, iters=5)
        print mtest, vtest    
    
    bkgf = SExtractorBackground(sigma_clip)
    bkg  = bkgf(nowdata[bmask])
    #mean, bkg, std = sigma_clipped_stats(nowdata[bmask], sigma=3.0, iters=5)
    #bkg   = np.nanmedian(nowdata[bmask])
    #print bkg, np.nanmean(nowdata[bmask]), np.nanstd(nowdata[bmask])
    
    apstar    = [CircularAperture([x[i],y[i]], r=r) for r in range(3,20)]
    nobkdata  = nowdata - bkg
    nobkdata[nobkdata < 0] = 0
    #anstar    = CircularAnnulus([x,y], r_in=5., r_out=10.)
    #apertures = [apstar, anstar]
    #print annulus
    phot_table = aperture_photometry(nobkdata, apstar)
    print '\n',phot_table,'\n', aperture_photometry(nowdata, apstar),'\n', bkg,'\n'
    #bkg_mean   = phot_table['aperture_sum_1'] / anstar.area()
    flux[i]    = phot_table['aperture_sum_2']*egain
    airmass[i] = airballoon.airmass(el*180.0/np.pi, 0) if el > 0 else np.nan#1.0/np.cos(np.pi/2.0 - el)
    dateflo[i] = float(loc.date)
    datetim.append(loc.date.datetime())
    #print '\n', loc.date, ephem.localtime(loc.date)
    #print '\n', float(loc.date), loc.date.datetime(), ephem.hours(loc.sidereal_time())
    #lst[i]     = ephem.hours(loc.sidereal_time())# * 12 / np.pi
    #print el, az, airmass[i]
#lst = np.array(lst)
#print datetim[-1]
#mask = (x > 1280) * (x < 0) * (y > 960) * (y < 0)

import matplotlib.pyplot as plt
mag = 2.5*np.log10(flux)
#mag[mask] = np.nan
fig, ax = plt.subplots(ncols=2, figsize=[10,4.5])
ma   = airmass < 4
scat = ax[0].scatter(airmass[ma], mag[ma], c=dateflo[ma], cmap='plasma')
cbar = plt.colorbar(scat, ax=ax[0])

print 'uno'

cbar.set_label('Time')
ax[0].set_xlabel('Airmass')
ax[0].set_ylabel(u'$2.5\log F$')
ax[0].set_title(star.name)

print 'dos'

#Star paths
if args.follow:
    stack = np.nanmedian(all_data[::5], axis=0)
    print 'tres'
    #fig, ax = plt.subplots()
    ax[1].matshow(stack, vmin=0, vmax=np.nanmedian(stack) + 2*np.nanstd(stack))
    ax[1].scatter(x[ma], y[ma], c=dateflo[ma], cmap='plasma')
    ax[1].set_title(folder_date)
    #aps = CircularAperture([x,y], r=5.)
    #aps.plot()
    #fig.savefig('%s_allsky_%s.png' % (args.starname, folder_date))
print 'aaa'
fig.tight_layout()
fig.savefig('%s_%s.png' % (star.name, folder_date))


#plt.show()
