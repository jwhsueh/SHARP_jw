import numpy as np
import matplotlib.pyplot as plt
import DistanceTool as distance
from astropy import convolution
import pyfits

HSTmock='../../data/illustris_1/mock_img/test_FID347491_SID0_ImageRGBLum.dat'

px_num=240
px_sc=0.04

bg=2.53445/1000. # counts/s
exp_t=4800 # 1 orbit=2400s
gain=7
inv_gain=1.0/gain

# zero pts
zpt=-21.1
wavelen=np.array([4297,5907,7995.937]) # A
fwhm=np.array([0.08,0.074,0.09]) # arcsec
photflam=2.508e-18 # inverse sensitivity, flux density one count per sec
perate=np.array([0,4096,2866]) #photoelectron rate
# approx 2 pixel

table=np.loadtxt(HSTmock)
F435W=table[:,0].reshape(px_num,px_num)
F606W=table[:,1].reshape(px_num,px_num)
F814W=table[:,2].reshape(px_num,px_num)

np.place(F435W,F435W==0,np.median(F435W))
np.place(F606W,F606W==0,np.median(F606W))
np.place(F814W,F814W==0,np.median(F814W))

## luminosity distance
class cosmopara:
	h = 0.704
	OM = 0.27

z=0.6
DL=distance.luminosity_distance(cosmopara,z)*1e6 # pc

m435=-2.5*np.log10(F435W)+5.*(np.log10(DL)-1)
m606=-2.5*np.log10(F606W)+5.*(np.log10(DL)-1)
m814=-2.5*np.log10(F814W)+5.*(np.log10(DL)-1)

## Vega magnitude
#m435=-2.5*np.log10(f435/Vf[0])

## PSF convolution

sig=2./2.3548
kernel=convolution.Gaussian2DKernel(sig)

#m814_cov=convolution.convolve(m814,kernel)
#np.place(m814,m814<=0,np.median(m814))
#print m814[:100]

## flux density

fv_814=3631*10**(-0.4*(m814)) # Jy
#np.place(fv_814,fv_814>100*np.median(fv_814),np.median(fv_814))
fl_814=fv_814*2.998e-5/wavelen[2]**2 # erg/cm2/s/A

ct_814=fl_814/photflam

# convolution
ct_814=convolution.convolve(ct_814,kernel)

# Poisson noise
pn_sig=np.sqrt(ct_814/inv_gain)

ct_814=ct_814*(1.0+np.random.normal(pn_sig))

#back to magnitude
mock_814fl=ct_814*photflam
mock_814fv=mock_814fl/2.998e-5*wavelen[2]**2

mock_814mag=-2.5*np.log10(mock_814fv/3631)


#hdulist=pyfits.PrimaryHDU(ct_814)
#hdulist.writeto('test2.fits')

plt.imshow(-mock_814mag)
plt.colorbar()
#plt.savefig('347491_F814W_poinoise.png')
plt.show()