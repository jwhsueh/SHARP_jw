from astropy.io import fits

hdulist=fits.open('./glafic/B1555_try6.fits')
model=hdulist[0].data

hdulist2=fits.open('./B1555/B1555_nirc2_n_Kp_6x6.fits')
image=hdulist2[0].data

res=image-model

hdu=fits.PrimaryHDU(res)
hdu.writeto('residue.fits')
