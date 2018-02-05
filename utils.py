import numpy as np

def altaz2xy(el, az):
    x0 = 1662.2182171489235
    y0 = 1231.5475985858097
    phi = 6.0238918377015276
    a0  = 812.65011866484531
    a1  = 12.110346726433246
    a2  = -41.790657061183573
    flipx = 0
    flipy = 1
    sizex = 3352
    sizey = 2532
    
    zen   = np.pi/2.0 - el
    theta = az + phi

    if (theta > 2*np.pi):
        theta = theta - 2*np.pi

    #theta = theta*np.pi/180.0

    r = a0*zen + a1*zen*zen + a2*zen*zen*zen

    x = x0 + r*np.cos(theta)
    y = r*np.sin(theta) + y0
    if flipx:
        x = sizex - x
    if flipy:
        y = sizey - y
        
    return x, y
