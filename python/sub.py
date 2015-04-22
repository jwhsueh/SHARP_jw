from astropy.io import fits

#hdulist=fits.open('./glafic/out_image.fits')
hdulist=fits.open('model_convolved.fits')
model=hdulist[0].data

#hdulist2=fits.open('./B1555/B1555_nirc2_n_Kp_6x6.fits')
hdulist2=fits.open('B1555_cutout_120x120pix.fits')

image=hdulist2[0].data

res=image-model

hdu=fits.PrimaryHDU(res)
hdu.writeto('Simona_residue.fits')
