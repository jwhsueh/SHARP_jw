from astropy.io import fits
from astropy.convolution import convolve
import numpy as np

#hdulist=fits.open('B1555_cimage_try13.fits')
#hdulist=fits.open('./B1555/B1555_nirc2_n_Kp_6x6.fits')
#hdulist=fits.open('/Users/jwhsueh/Downloads/bgsub_17.fits')
hdulist=fits.open('./B0712_17_gal.fits')
image=hdulist[0].data

box = np.zeros((5,5))
box.fill(1./25.)
box[2,2]=0


#cutout zone (1" view)
#s=s[255:356,249:350] #for SHARP image

#s=s[50:250,50:250]
#s1=image[575:675,79:179] #for psf image
#s=image[149:349,617:817] # for galaxy


## ---convolution

cimage=convolve(image,box)
med = np.median(image)
print med

diff = np.abs(image)
flag = [diff>10.*med]
#print image[flag]
#print cimage[flag]

image[flag] = cimage[flag]

#hdu1=fits.PrimaryHDU(cimage)

#hdu1=fits.PrimaryHDU(s1)
hdu=fits.PrimaryHDU(image)
#hdu.writeto('B1555_model_cutout.fits')
#hdu.writeto('B1555_Kp_cutout.fits')
#hdu1.writeto('B0712_01_c_rdx2.fits',clobber = True)
hdu.writeto('B0712_17_gal_rdx.fits',clobber = True)
#hdu1.writeto('B0712_17_psf.fits',clobber = True)
#hdu.writeto('B0712_17_gal.fits',clobber = True)

