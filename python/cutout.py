from astropy.io import fits

#hdulist=fits.open('B1555_cimage_try13.fits')
#hdulist=fits.open('./B1555/B1555_nirc2_n_Kp_6x6.fits')
hdulist=fits.open('SIE2_per_image.fits')
s=hdulist[0].data

#cutout zone (1" view)
#s=s[255:356,249:350] #for SHARP image

s=s[205:406,199:400] #for SHARP image

hdu=fits.PrimaryHDU(s)
#hdu.writeto('B1555_model_cutout.fits')
#hdu.writeto('B1555_Kp_cutout.fits')
hdu.writeto('B1555_glafic_cutout.fits')
