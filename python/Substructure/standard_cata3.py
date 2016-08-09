import numpy as np
import h5py
import DistanceTool as distance
import snapshot

import pandas as pd

basePath = '../../data/illustris_1'
redshift = 0.2

ssNumber = '120'
print ssNumber


## maybe we don't need this part
'''
## ----------- Dandan's morphology pick --------------- ##
# read-in galaxy IDs
DanM_fx,DanM_fy,DanM_fz = ['subfindID','R_Ex'],['subfindID','R_Ey'],['subfindID','R_Ez']
DanM_x = pd.read_csv('JenWei_table_disk_in_x.dat',sep = '\t',names = DanM_fx,comment = '#')
DanM_y = pd.read_csv('JenWei_table_disk_in_y.dat',sep = '\t',names = DanM_fy,comment = '#')
DanM_z = pd.read_csv('JenWei_table_disk_in_z.dat',sep = '\t',names = DanM_fz,comment = '#')

DanM = pd.merge(DanM_x,DanM_y,how = 'inner',on = 'subfindID')
DanM = pd.merge(DanM,DanM_z,how = 'inner',on = 'subfindID')

DanM = DanM.sort(['subfindID'])
DanM['subfindID'] = DanM['subfindID'].astype(int)

DanM = DanM['subfindID']
'''
## --------- kinematics ------------ ##

catalog = basePath+'/kinematics_'+str(ssNumber)+'_test.dat'
kine_field = ['subfindID','disk_star_f','bulge_star_f']
kinematics = pd.read_csv(catalog,sep = '\s+',names = kine_field,comment = '#')

## -------- group catalog properties [general] ----- ##

catalog2 = basePath+'/Galaxy_'+str(ssNumber)+'_test.dat'
group_field = ['subfindID','pos_x','pos_y','pos_z','mass','stellar_mass_all','velDisp','R_halfmass','relaxation','r_mag','B_mag','V_mag']
group = pd.read_csv(catalog2,sep = '\s+',names = group_field,comment = '#')

# standard DataFrame
standard = pd.merge(group,kinematics,how = 'inner',on = 'subfindID')
print np.array(group.subfindID)[-1],np.array(standard.subfindID)[-1]

'''
## -------- inclination angle --------- ##

catalog = basePath+'/inclination_'+ssNumber+'_sig.dat'

IDs = np.loadtxt(catalog,dtype = 'int',unpack=True, usecols=[0])
theta_x = np.loadtxt(catalog,dtype = 'float',unpack=True, usecols=[1])
theta_y = np.loadtxt(catalog,dtype = 'float',unpack=True, usecols=[2])
theta_z = np.loadtxt(catalog,dtype = 'float',unpack=True, usecols=[3])

inc_x = pd.DataFrame({'subfindID':IDs,'theta':theta_x})
inc_y = pd.DataFrame({'subfindID':IDs,'theta':theta_y})
inc_z = pd.DataFrame({'subfindID':IDs,'theta':theta_z})
'''
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

#Lens_cata = pd.merge(Lens_x,Lens_y,how='inner',on = 'subfindID')
#Lens_cata = pd.merge(Lens_cata,Lens_z,how='inner',on = 'subfindID')

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

Phot_x = pd.DataFrame(Phot_x,columns = ['subfindID','stellar_mass','Sersic','SteR_halfmass','SB_Exp','morphology'])
Phot_y = pd.DataFrame(Phot_y,columns = ['subfindID','stellar_mass','Sersic','SteR_halfmass','SB_Exp','morphology'])
Phot_z = pd.DataFrame(Phot_z,columns = ['subfindID','stellar_mass','Sersic','SteR_halfmass','SB_Exp','morphology'])

#Phot_cata = pd.merge(Phot_x,Phot_y,how='inner',on = 'subfindID')
#Phot_cata = pd.merge(Phot_cata,Phot_z,how='inner',on = 'subfindID')

#print Phot_cata

## ----- FULL gas catalog (Dandan) ----- ##

catalog3x = basePath+'/Dandan_Gas'+str(ssNumber)+'_x.dat'
catalog3y = basePath+'/Dandan_Gas'+str(ssNumber)+'_y.dat'
catalog3z = basePath+'/Dandan_Gas'+str(ssNumber)+'_z.dat'

Gas_fx = ['subfindID','fgas','fcgs']
Gas_fy = ['subfindID','fgas','fcgs']
Gas_fz = ['subfindID','fgas','fcgs']

Gas_x = pd.read_csv(catalog3x,sep = '\s+',names = Gas_fx,comment = '#')
Gas_y = pd.read_csv(catalog3y,sep = '\s+',names = Gas_fy,comment = '#')
Gas_z = pd.read_csv(catalog3z,sep = '\s+',names = Gas_fz,comment = '#')



## -------- put together to larger cataloga (x,y,z) ------- ##

standard_x = pd.merge(standard,Lens_x,how = 'inner',on = 'subfindID')
#standard_x = pd.merge(standard_x,inc_x,how = 'inner',on = 'subfindID')
standard_x = pd.merge(standard_x,Phot_x,how = 'inner',on = 'subfindID')
standard_x = pd.merge(standard_x,Gas_x,how = 'inner',on = 'subfindID')

standard_y = pd.merge(standard,Lens_y,how = 'inner',on = 'subfindID')
#standard_y = pd.merge(standard_y,inc_y,how = 'inner',on = 'subfindID')
standard_y = pd.merge(standard_y,Phot_y,how = 'inner',on = 'subfindID')
standard_y = pd.merge(standard_y,Gas_y,how = 'inner',on = 'subfindID')

standard_z = pd.merge(standard,Lens_z,how = 'inner',on = 'subfindID')
#standard_z = pd.merge(standard_z,inc_z,how = 'inner',on = 'subfindID')
standard_z = pd.merge(standard_z,Phot_z,how = 'inner',on = 'subfindID')
standard_z = pd.merge(standard_z,Gas_z,how = 'inner',on = 'subfindID')


#standard = pd.merge(standard,Lens_cata,how = 'inner',on = 'subfindID')
#standard = pd.merge(standard,Phot_cata,how = 'inner',on = 'subfindID')

#print standard

standard_x.to_csv(basePath+'/test_'+str(ssNumber)+'_x.dat',sep = '\t',index = False)
standard_y.to_csv(basePath+'/test_'+str(ssNumber)+'_y.dat',sep = '\t',index = False)
standard_z.to_csv(basePath+'/test_'+str(ssNumber)+'_z.dat',sep = '\t',index = False)
