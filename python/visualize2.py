import numpy as np
import matplotlib.pyplot as plt
#import pyfits
from astropy.io import fits

re=100
L=1.5 #arcsec
pix=L/re #pixel size

x=np.linspace(-L/2.0,L/2.0,re)
y=np.linspace(-L/2.0,L/2.0,re)

xi,yi=np.meshgrid(x,y)

#parameters
##SIE
b=0.2
e=0.56
q=1-e
esp=np.sqrt((1.0-q)/(1.0+q))
c=[-1.651000e-01, -2.403000e-01] #center of fig

PA=-3.0
st=90.0 #set east to north
angle=np.radians((PA)+st)

##expdisk
k0=0.5
rd=0.105
e2=0.85
q2=1.0-e2
c2=[-1.7e-01, -2.0e-01]

PA2=-10.0
angle2=np.radians((PA2)+st)

#translation for expdisk
x2=x+(c[0]-c2[0])
y2=y+(c[1]-c2[1])

x2i,y2i=np.meshgrid(x2,y2)

#inverse rotate

vx=np.cos(angle)*xi+np.sin(angle)*yi
vy=-1.0*np.sin(angle)*xi+np.cos(angle)*yi

vx2=np.cos(angle2)*x2i+np.sin(angle2)*y2i
vy2=-1.0*np.sin(angle2)*x2i+np.cos(angle2)*y2i

kappa=np.zeros((re,re))
kappa2=np.zeros((re,re))
for i in range(re):
	for j in range(re):
		kappa[j,i]=0.5*b/np.sqrt((1.0-esp)*vx[j,i]**2+(1.0+esp)*vy[j,i]**2)
		kappa2[j,i]=k0*np.exp(-1.0*np.sqrt(vx2[j,i]**2+vy2[j,i]**2/q2**2)/rd)/q


##fits

hdu=fits.PrimaryHDU(kappa+kappa2)
hdr=hdu.header

#CD matrix

off=[0.0-c[0],0.0-c[1]]
x0,y0=np.interp(off[0],x,np.arange(re)),np.interp(off[1],y,np.arange(re)) #interploration
a0,d0=239.3167,37.36

c11,c22=pix,pix
c12,c21=0.0,0.0

#hdr.set(('CD1_1','CD1_2','CD2_1','CD2_2','CRPIX1','CRPIX2','CRVAL1','CRVAL2'))

hdr['CD1_1'],hdr['CD1_2']=c11,c12
hdr['CD2_1'],hdr['CD2_2']=c21,c22

hdr['CRPIX1'],hdr['CRPIX2']=x0,y0
hdr['CRVAL1'],hdr['CRVAL2']=a0,d0

hdr['CTYPE1'],hdr['CTYPE2']='RA---TAN','DEC--TAN'

hdu.writeto('mass_profile.fits')

##plot
#fig1=plt.imshow(np.log(kappa+kappa2),origin='bottom')
#plt.clim(0.,50)
#plt.colorbar()
#plt.show()
