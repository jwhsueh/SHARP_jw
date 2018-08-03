import numpy as np
import snapshot
import h5py
import groupcat
import DistanceTool as distance
import matplotlib.pyplot as plt
from astropy.io import fits
from matplotlib.colors import LogNorm


class cosmopara:
	h = 0.704
	OM = 0.27

snapshotPath = '/Volumes/narsil_1/jwhsueh/illustris_1'
snapNum = '089'
ssNum = 89
galaxy_list = np.loadtxt(snapshotPath+'/snap89_reid.txt')
galaxy_list = galaxy_list.astype(int)

for GalaxyID in galaxy_list:

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
	#SFT = np.sort(SFT)
	#print SFT[:100]

	## change ref point to subhalo CM
	st_x,st_y,st_z = st_x-CM_x,st_y-CM_y,st_z-CM_z	# ckpc/h
	st_x,st_y,st_z = st_x*a,st_y*a,st_z*a # kpc/h

	#mask = SFT>0.60
	mask = SFT>np.max(SFT)*0.7
	print np.max(SFT)*0.7

	st_x,st_y,st_z = st_x[mask],st_y[mask],st_z[mask]

	#plt.scatter(st_z,-st_y,marker='.',s=1) #proj x
	plt.hist2d(st_z,-st_y,bins=100,norm=LogNorm()) #proj x
	plt.savefig('../../data/illustris_1/snapshot_'+str(ssNum)+'_fig/subfind_'+str(GalaxyID)+'_x.png')
	plt.clf()
	plt.hist2d(st_x,st_y,bins=100,norm=LogNorm())  #proj z
	plt.savefig('../../data/illustris_1/snapshot_'+str(ssNum)+'_fig/subfind_'+str(GalaxyID)+'_z.png')
	plt.clf()
	plt.hist2d(st_x,st_z,bins=100,norm=LogNorm())  #proj y
	plt.savefig('../../data/illustris_1/snapshot_'+str(ssNum)+'_fig/subfind_'+str(GalaxyID)+'_y.png')
	plt.clf()
	
	#plt.show()

'''
## critical curve part
hdu = fits.open('/Volumes/sting_1/snap99_result/snap99_'+subfindID+'/critical_'+subfindID+'_p3_64NN.fits')
data = hdu[0].data

crit_y,crit_x = np.where(data==2)
crit_x,crit_y = crit_x-512,crit_y-512
pix_size = 0.00206 # arcsec
pix_size = pix_size*6.730/cosmopara.h#*1.2 # kpc/h
plt.scatter(crit_x*pix_size,crit_y*pix_size,color='r',marker='.',s=1)
plt.xlim(-12,12)
plt.ylim(-12,12)
plt.xlabel('kpc/h',fontsize=14)
#plt.ylabel('kpc/h',fontsize=14)
plt.text(-11,10.5,'subhaloID:'+subfindID,fontsize=16)
#plt.title('subfind '+subfindID+' proj_3')
plt.gca().set_aspect('equal')
plt.show()
#plt.savefig('/Users/jwhsueh/Documents/SHARP_jw/data/illustris_1/mock_img/'+subfindID+'_ell.png',bbox_inches='tight')
'''