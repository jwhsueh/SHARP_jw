from astropy.io import fits
from astropy.convolution import convolve

hdulist=fits.open('B1422_PSF.fits')
psf=hdulist[0].data

psf=psf[1:200,1:200]

#hdulist2=fits.open('./glafic/B1555_try11.fits')
hdulist2=fits.open('B1555_glafic_cutout.fits')
image=hdulist2[0].data

cimage=convolve(image,psf)

hdu=fits.PrimaryHDU(cimage)
hdu.writeto('B1555_cimage.fits')
