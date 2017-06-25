from astropy.io import fits
from astropy import convolution
from astropy import wcs
import matplotlib.pyplot as plt
import numpy as np
from scipy import misc

hdulist=fits.open('/Volumes/sting_1/snap99_140304/image_140304_p1_64src_10.fits')
image=hdulist[0].data
w=wcs.WCS(hdulist[0].header)

real_size = 0.0017491337 # deg
real_size = real_size*3600 # arcsec
src = np.array([[-0.000344586,0.000394406]]) # rad
src = src/real_size*3600*2/1.1
print src


x = np.linspace(-1*real_size/2.,real_size/2.,256)

#psf=psf[1:200,1:200]
kernel=convolution.Gaussian2DKernel(3)
image=convolution.convolve(image,kernel)

hdulist2=fits.open('/Volumes/sting_1/snap99_140304/critical_140304_p1_64NN.fits')
array=hdulist2[0].data

xx = np.linspace(-1*real_size/2/1.2,real_size/2/1.2,1024)

plt.contour(xx,xx,array,colors='r')
#plt.show()

hdulist3=fits.open('/Volumes/sting_1/snap99_140304/caustics_140304_p1_64NN.fits')
array2=hdulist3[0].data

plt.contour(xx,xx,array2,colors='b')

##

plt.contour(x,x,image,colors='k')
plt.scatter(-1*src[0][0],-1*src[0][1],marker='o',edgecolor='k',color='g',s=70)
plt.xlim(-2.5,2.5)
plt.ylim(-2.5,2.5)
plt.gca().set_aspect('equal')
plt.xlabel(r'$\Delta \alpha (arcsec)$',fontsize=16)
plt.ylabel('$\Delta \delta (arcsec)$',fontsize=16)
#plt.show()
plt.savefig('../data/glamer/140304_p1_64src_10_lens.png',bbox_inches='tight')
