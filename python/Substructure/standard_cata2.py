import numpy as np
import h5py
import DistanceTool as distance
import snapshot

import pandas as pd

basePath = '../../data/illustris_1'
ssNumber = 99
Snapshot_num = 'Snapshot_'+str(ssNumber)

## get header file
headerPath = '/Volumes/narsil_1/jwhsueh/illustris_1'

f = h5py.File(snapshot.snapPath(headerPath,ssNumber),'r')
header = dict( f['Header'].attrs.items() )

redshift = header['Redshift']


## maybe we don't need this part
'''
## ----------- Dandan's morphology pick --------------- ##
# read-in galaxy IDs
DanM_fx,DanM_fy,DanM_fz = ['subfindID','R_Ex'],['subfindID','R_Ey'],['subfindID','R_Ez']
DanM_x = pd.read_csv('JenWei_table_disk_in_x.dat',sep = '\t',names = DanM_fx,comment = '#')
DanM_y = pd.read_csv('JenWei_table_disk_in_y.dat',sep = '\t',names = DanM_fy,comment = '#')
DanM_z = pd.read_csv('JenWei_table_disk_in_z.dat',sep = '\t',names = DanM_fz,comment = '#')

DanM = pd.merge(DanM_x,DanM_y,how = 'outer',on = 'subfindID')
DanM = pd.merge(DanM,DanM_z,how = 'outer',on = 'subfindID')

DanM = DanM.sort(['subfindID'])
DanM['subfindID'] = DanM['subfindID'].astype(int)

DanM = DanM['subfindID']
'''
## --------- kinematics ------------ ##

catalog = basePath+'/kinematics_'+str(ssNumber)+'_sig.dat'
kine_field = ['subfindID','disk_star_f','bulge_star_f']
kinematics = pd.read_csv(catalog,sep = '\s+',names = kine_field,comment = '#')

## -------- group catalog properties ----- ##

catalog2 = basePath+'/Galaxy_'+str(ssNumber)+'_sig.dat'
group_field = ['subfindID','pos_x','pos_y','pos_z','mass','stellar_mass','velDisp','R_halfmass','relaxation','r_mag']
group = pd.read_csv(catalog2,sep = '\s+',names = group_field,comment = '#')

# standard DataFrame
standard = pd.merge(group,kinematics,how = 'inner',on = 'subfindID')

## ----- FULL lens catalog (Dandan) ----- ##

catalog3x = basePath+'/Dandan_Lens'+str(ssNumber)+'_x.dat'
catalog3y = basePath+'/Dandan_Lens'+str(ssNumber)+'_y.dat'
catalog3z = basePath+'/Dandan_Lens'+str(ssNumber)+'_z.dat'

Lens_f = ['subfindID','mass','2dmass_R_E','R_E','DMfrac_R_E']
#Lens_fy = ['subfindID','mass','2dmass_R_Ey','R_Ey','DMfrac_R_Ey']
#Lens_fz = ['subfindID','mass','2dmass_R_Ez','R_Ez','DMfrac_R_Ez']

Lens_x = pd.read_csv(catalog3x,sep = '\s+',names = Lens_f,comment = '#')
Lens_y = pd.read_csv(catalog3y,sep = '\s+',names = Lens_f,comment = '#')
Lens_z = pd.read_csv(catalog3z,sep = '\s+',names = Lens_f,comment = '#')

Lens_x = pd.DataFrame(Lens_x,columns = ['subfindID','2dmass_R_E','R_E','DMfrac_R_E'])
Lens_y = pd.DataFrame(Lens_y,columns = ['subfindID','2dmass_R_E','R_E','DMfrac_R_E'])
Lens_z = pd.DataFrame(Lens_z,columns = ['subfindID','2dmass_R_E','R_E','DMfrac_R_E'])

#Lens_cata = pd.merge(Lens_x,Lens_y,how='outer',on = 'subfindID')
#Lens_cata = pd.merge(Lens_cata,Lens_z,how='outer',on = 'subfindID')

#print Lens_cata

## ----- FULL photometry catalog (Dandan) ----- ##

catalog3x = basePath+'/Dandan_Phot'+str(ssNumber)+'_x.dat'
catalog3y = basePath+'/Dandan_Phot'+str(ssNumber)+'_y.dat'
catalog3z = basePath+'/Dandan_Phot'+str(ssNumber)+'_z.dat'

Phot_fx = ['subfindID','mass','stellar_mass','Sersic','SteR_halfmass','SB_Exp','morphology']
Phot_fy = ['subfindID','mass','stellar_mass','Sersic','SteR_halfmass','SB_Exp','morphology']
Phot_fz = ['subfindID','mass','stellar_mass','Sersic','SteR_halfmass','SB_Exp','morphology']

Phot_x = pd.read_csv(catalog3x,sep = '\s+',names = Phot_fx,comment = '#')
Phot_y = pd.read_csv(catalog3y,sep = '\s+',names = Phot_fy,comment = '#')
Phot_z = pd.read_csv(catalog3z,sep = '\s+',names = Phot_fz,comment = '#')

Phot_x = pd.DataFrame(Phot_x,columns = ['subfindID','Sersic','SteR_halfmass','SB_Exp','morphology'])
Phot_y = pd.DataFrame(Phot_y,columns = ['subfindID','Sersic','SteR_halfmass','SB_Exp','morphology'])
Phot_z = pd.DataFrame(Phot_z,columns = ['subfindID','Sersic','SteR_halfmass','SB_Exp','morphology'])

#Phot_cata = pd.merge(Phot_x,Phot_y,how='outer',on = 'subfindID')
#Phot_cata = pd.merge(Phot_cata,Phot_z,how='outer',on = 'subfindID')

#print Phot_cata

## -------- put together to larger cataloga (x,y,z) ------- ##

standard_x = pd.merge(standard,Lens_x,how = 'inner',on = 'subfindID')
standard_x = pd.merge(standard_x,Phot_x,how = 'inner',on = 'subfindID')

standard_y = pd.merge(standard,Lens_y,how = 'inner',on = 'subfindID')
standard_y = pd.merge(standard_y,Phot_y,how = 'inner',on = 'subfindID')

standard_z = pd.merge(standard,Lens_z,how = 'inner',on = 'subfindID')
standard_z = pd.merge(standard_z,Phot_z,how = 'inner',on = 'subfindID')


#standard = pd.merge(standard,Lens_cata,how = 'inner',on = 'subfindID')
#standard = pd.merge(standard,Phot_cata,how = 'inner',on = 'subfindID')

#print standard

standard_x.to_csv(basePath+'/full_'+str(ssNumber)+'_x.dat',sep = '\t',index = False)
standard_y.to_csv(basePath+'/full_'+str(ssNumber)+'_y.dat',sep = '\t',index = False)
standard_z.to_csv(basePath+'/full_'+str(ssNumber)+'_z.dat',sep = '\t',index = False)
