import numpy as np
import snapshot
import h5py
import groupcat

snapshotPath = '/Volumes/narsil_1/jwhsueh/illustris_1'
snapNum = '099'
ssNum = 99
GalaxyID = 139175
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
subhalo_dm = snapshot.loadSubhalo(snapshotPath,snapNum,GalaxyID,'dm')
subhalo_gas = snapshot.loadSubhalo(snapshotPath,snapNum,GalaxyID,'gas')

st_ms = subhalo_st['Masses']*1e10/0.704
gas_ms = subhalo_gas['Masses']*1e10/0.704

coord = subhalo_st['Coordinates']
st_x,st_y,st_z = coord[:,0],coord[:,1],coord[:,2]