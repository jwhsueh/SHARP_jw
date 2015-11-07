## Now it's only for change the scale size

import pyfits
import numpy as np
from astropy.convolution import convolve

ratio=5 # rescale ratio between pixel sizes
root='/Users/jwhsueh/Documents/SHARP_jw/models/' # root path for fits files
system='B0712' # name of system
file='B0712_nirc2_n_Kp_6x6.fits'
kfile='nirc2nic_GaussKernel.fits' # kernel is for original image

ori=pyfits.open(root+system+'/'+file)
orimg=ori[0].data


def resamp(img):
    #ori=pyfits.open('/Users/jwhsueh/Documents/SHARP_jw/models/B1555/B1555_nirc2_4x4.fits')
    #ori=pyfits.open(root+system+'/'+file)
    #img=ori[0].data

    size=img.shape[0]

    new=np.zeros((size/ratio,size/ratio))
    ace=np.zeros((ratio,ratio))

    for i in range(size/ratio):
        for j in range(size/ratio):
            ace=img[ratio*i:(ratio*i+ratio),ratio*j:(ratio*j+ratio)]
            index=np.sum(ace)/ratio**2

            new[i,j]=index

    return new


#print new.shape

    #hdu=pyfits.PrimaryHDU(new)
    #hdu.writeto('/Users/jwhsueh/Documents/SHARP_jw/models/B1555/test.fits')

def conv(img):
    kernel_file=pyfits.open(root+kfile)
    kernel=kernel_file[0].data

    new=convolve(img,kernel)

    return new


new=conv(orimg)
new=resamp(new)

hdu=pyfits.PrimaryHDU(new)
hdu.writeto(root+system+'/test.fits')


