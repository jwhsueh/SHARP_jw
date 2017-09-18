import numpy as np
import snapshot
import h5py
import groupcat
import DistanceTool as distance
import matplotlib.pyplot as plt
from astropy.io import fits


class cosmopara:
	h = 0.704
	OM = 0.27

snapshotPath = '/Volumes/narsil_1/jwhsueh/illustris_1'
snapNum = '099'
ssNum = 99
GalaxyID = 154402
subfindID = str(GalaxyID)


f = h5py.File(snapshot.snapPath(snapshotPath,ssNum),'r')
header = dict( f['Header'].attrs.items() )
redshift = header['Redshift']

a = 1.0/(1.0+redshift) # scale factor

# get peculiar velocity of group
SubPos = groupcat.loadSubhalos(snapshotPath,snapNum, fields = ['SubhaloPos'])[GalaxyID] # ckpc/h
CM_x,CM_y,CM_z = SubPos[0],SubPos[1],SubPos[2]
print CM_x

subhalo_st = snapshot.loadSubhalo(snapshotPath,snapNum,GalaxyID,'stars')

coord = subhalo_st['Coordinates']
st_x,st_y,st_z = coord[:,0],coord[:,1],coord[:,2]

SFT = subhalo_st['GFM_StellarFormationTime'] # in scale factor

## change ref point to subhalo CM
st_x,st_y,st_z = st_x-CM_x,st_y-CM_y,st_z-CM_z	# ckpc/h
st_x,st_y,st_z = st_x*a,st_y*a,st_z*a # kpc/h

#mask = SFT>0.61
mask = SFT>0.5

st_x,st_y,st_z = st_x[mask],st_y[mask],st_z[mask]

plt.scatter(st_x,st_y,marker='.',s=1)
#plt.show()


## critical curve part
hdu = fits.open('/Volumes/sting_1/snap99_result/snap99_'+subfindID+'/critical_'+subfindID+'_p3_64NN.fits')
data = hdu[0].data

crit_y,crit_x = np.where(data==2)
crit_x,crit_y = crit_x-512,crit_y-512
pix_size = 0.00206 # arcsec
pix_size = pix_size*6.730/cosmopara.h*1.2 # kpc/h
plt.scatter(crit_x*pix_size,crit_y*pix_size,color='r',s=1)
plt.xlim(-20,20)
plt.ylim(-20,20)
plt.xlabel('kpc/h')
plt.ylabel('kpc/h')
plt.title('subfind '+subfindID+' proj_3')
plt.gca().set_aspect('equal')
#plt.show()
plt.savefig(subfindID+'_ell.png')