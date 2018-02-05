from photutils import DAOStarFinder, CircularAperture
from astropy.stats import sigma_clipped_stats, sigma_clip
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import sys

data = fits.getdata(sys.argv[1])
#hdr  = fits.getheader(sys.argv[1], 1)
#print hdr
#egain = hdr['EGAIN']
#ped   = hdr['PEDESTAL']
#data  = (data + ped) * egain
data = sigma_clip(data, sigma=3, iters=3)
print data.max(), data.min()

mean, median, std = sigma_clipped_stats(data, sigma=3.0, iters=5)
print mean,median,std

daofind = DAOStarFinder(fwhm=3.0, threshold=3*std)
sources = daofind(data)

x, y = sources['xcentroid'], sources['ycentroid']
apertures  = CircularAperture([x,y], r=4.)

fig, ax = plt.subplots()
ax.imshow(data, cmap='Greys_r')
apertures.plot(color='lime', lw=1.5, alpha=0.5)

plt.show()
