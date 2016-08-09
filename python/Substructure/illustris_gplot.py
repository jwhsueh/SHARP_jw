import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

basePath = '../../data/illustris_1/'

ssNumber = ['085','099','120'] # snapshotnumber
z = [1.0,0.6,0.2]

## -------- group catalog properties [general] ----- ##

catalog = basePath+'/Galaxy_'+ssNumber[0]+'_test.dat'
group_field = ['subfindID','pos_x','pos_y','pos_z','mass','stellar_mass','velDisp','R_halfmass','relaxation','r_mag','B_mag','V_mag']
group0 = pd.read_csv(catalog,sep = '\s+',names = group_field,comment = '#')

catalog = basePath+'/Galaxy_'+ssNumber[1]+'_test.dat'
group1 = pd.read_csv(catalog,sep = '\s+',names = group_field,comment = '#')

catalog = basePath+'/Galaxy_'+ssNumber[2]+'_test.dat'
group2 = pd.read_csv(catalog,sep = '\s+',names = group_field,comment = '#')

## stellar_mass cut

z0 = group0.loc[(group0.stellar_mass>1e9),:]
z1 = group1.loc[(group1.stellar_mass>1e9),:]
z2 = group2.loc[(group2.stellar_mass>1e9),:]

z0['log_stellar_mass'] = np.log10(z0.stellar_mass)
z1['log_stellar_mass'] = np.log10(z1.stellar_mass)
z2['log_stellar_mass'] = np.log10(z2.stellar_mass)


## scatter plot

plt.scatter(z0.velDisp,z0.log_stellar_mass,c = 'b',label = 'z='+str(z[0]),marker = 'x')
plt.scatter(z1.velDisp,z1.log_stellar_mass,c='g',label = 'z='+str(z[1]),marker = 'x')
plt.scatter(z2.velDisp,z2.log_stellar_mass,c = 'r',label = 'z='+str(z[2]),marker = 'x')
plt.title('All galaxies')
plt.show()

'''
se = np.linspace(0,370,40)
dot = []
for i in range(se.size-1):
	dot.append((se[i]+se[i+1])/2.)

group0_vel = np.histogram(group0.velDisp,bins = se)[0].astype(float)
group1_vel = np.histogram(group1.velDisp,bins = se)[0].astype(float)
group2_vel = np.histogram(group2.velDisp,bins = se)[0].astype(float)

## plot

plt.plot(dot,group0_vel,color='k',drawstyle='steps',label = 'z='+str(z[0]))
plt.plot(dot,group1_vel,color='b',drawstyle='steps',label = 'z='+str(z[1]))
plt.plot(dot,group2_vel,color='r',drawstyle='steps',label = 'z='+str(z[2]))
plt.legend()
plt.xlabel('velocity dispersion [km/s]')
plt.ylabel('counts')
plt.title('All galaxies')
#plt.show()
plt.savefig(basePath+'redshift_evo/velDisp.png')
'''