import pandas as pd
import seaborn as sns
import numpy as np

import matplotlib.pyplot as plt
import DistanceTool as distance

basePath = '../../data/illustris_1/redshift_evo/'
ssNumber = ['085','099','120']
redshift = [1.0,0.6,0.2]

file0 = basePath+'snapshot'+ssNumber[0]+'st_hist_30.dat'
file1 = basePath+'snapshot'+ssNumber[1]+'st_hist_30.dat'
file2 = basePath+'snapshot'+ssNumber[2]+'st_hist_30.dat'

file0_all = basePath+'snapshot'+ssNumber[0]+'st_hist_all.dat'
file1_all = basePath+'snapshot'+ssNumber[1]+'st_hist_all.dat'
file2_all = basePath+'snapshot'+ssNumber[2]+'st_hist_all.dat'
'''
## scatter plot
z0,z1,z2 = pd.read_csv(file0,sep = '\t'),pd.read_csv(file1,sep = '\t'),pd.read_csv(file2,sep = '\t')
#print z0

plt.scatter(z0.V_mag,z0.color,c = 'b',label = 'z='+str(redshift[0]),marker = '.')
plt.scatter(z1.V_mag,z1.color,c='g',label = 'z='+str(redshift[1]),marker = 'x')
plt.scatter(z2.V_mag,z2.color,c = 'r',label = 'z='+str(redshift[2]))

plt.title('Galaxies satisfy all criteria')
plt.legend(loc = 1)
plt.xlabel('V_mag')
plt.ylabel('B-V')
plt.xlim(15,35)
plt.ylim(-0.2,1.)
#plt.show()
plt.savefig('./illustris_cross/zEvo_color_three.png',bbox_inches='tight')

'''
## hist

z0,z1,z2 = pd.read_csv(file0,sep = '\t'),pd.read_csv(file1,sep = '\t'),pd.read_csv(file2,sep = '\t')
z0_all,z1_all,z2_all = pd.read_csv(file0_all,sep = '\t'),pd.read_csv(file1_all,sep = '\t'),pd.read_csv(file2_all,sep = '\t')

z0['mor_hist'] = z0.morphology/z0.All
z1['mor_hist'] = z1.morphology/z1.All
z2['mor_hist'] = z2.morphology/z2.All

z0['disk_hist'] = z0.disk_star/z0.All
z1['disk_hist'] = z1.disk_star/z1.All
z2['disk_hist'] = z2.disk_star/z2.All

z0['bulge_hist'] = z0.bulge_star/z0.All
z1['bulge_hist'] = z1.bulge_star/z1.All
z2['bulge_hist'] = z2.bulge_star/z2.All

#print z0

plt.plot(z0.log_stellar_mass,z0.disk_hist,color='k',label = 'z='+str(redshift[0]),drawstyle='steps',linestyle = 'dashdot')
plt.plot(z0.log_stellar_mass,z1.disk_hist,'b',label = 'z='+str(redshift[1]),drawstyle='steps',linestyle = 'dashed')
plt.plot(z0.log_stellar_mass,z2.disk_hist,'r',label = 'z='+str(redshift[2]),drawstyle='steps')

plt.xlabel('log(stellar mass)')
plt.ylabel('Disk System Fraction')
#plt.title('morphology')
#plt.title('Bulge star < 60%')
plt.title('Thin disk star > 40%')
plt.legend(scatterpoints=1,loc =1)
#plt.legend(loc = 1)
plt.ylim(0,0.8)
plt.xlim(10,12)
#plt.show()
plt.savefig('./illustris_cross/zEvo_st30_hist_df_2.png',bbox_inches='tight')
