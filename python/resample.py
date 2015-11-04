## Now it's only for change the scale size

import pyfits
import numpy as np

ratio=5 # rescale ratio between pixel sizes

ori=pyfits.open('/Users/jwhsueh/Documents/SHARP_jw/models/B1555/B1555_nirc2_4x4.fits')
img=ori[0].data

size=img.shape[0]

new=np.zeros((size/ratio,size/ratio))
ace=np.zeros((ratio,ratio))

for i in range(size/ratio):
    for j in range(size/ratio):
        ace=img[ratio*i:(ratio*i+ratio),ratio*j:(ratio*j+ratio)]
        index=np.sum(ace)/ratio**2

        new[i,j]=index


#print new.shape

hdu=pyfits.PrimaryHDU(new)
hdu.writeto('/Users/jwhsueh/Documents/SHARP_jw/models/B1555/test.fits')