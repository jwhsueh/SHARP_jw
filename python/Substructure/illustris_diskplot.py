import pandas as pd
import seaborn as sns
import numpy as np

import matplotlib.pyplot as plt
import DistanceTool as distance

basePath = '../../data/illustris_1/'
ssNumber = '120'
redshift = 0.2
print ssNumber

catalog_x = basePath+'test_'+str(ssNumber)+'_x.dat'
catalog_y = basePath+'test_'+str(ssNumber)+'_y.dat'
catalog_z = basePath+'test_'+str(ssNumber)+'_z.dat'
table_x,table_y,table_z = pd.read_csv(catalog_x,sep = '\t'),pd.read_csv(catalog_y,sep = '\t'),pd.read_csv(catalog_z,sep = '\t')

## to a long table
table = table_x.append(table_y,ignore_index = True)
all_table = table.append(table_z,ignore_index = True)

all_table = all_table.loc[(all_table.r_mag <100)&(all_table.B_mag <100)&(all_table.V_mag <100),:]

#print all_table

## stellar mass lower limit
all_table = all_table.loc[(all_table.stellar_mass>1e9),:]

## add columns

all_table.stellar_mass = all_table.stellar_mass/0.71
all_table['log_stellar_mass']= np.log10(all_table.stellar_mass)
all_table['color'] = all_table.B_mag - all_table.V_mag
all_table['mor_mix'] = 0

## total flux cut

#all_table = all_table.loc[(all_table.r_mag<25.),:]


## relaxation
table = all_table.loc[(all_table.relaxation == 0),:] # table is for disk systems
#table = all_table # disk_table

## edge-on pick
#table = table.loc[(table.theta > 80.) & (table.theta < 100.),:]



#print table.log_stellar_mass

## add apparent surface brightness
''' cosmopara '''
class cosmopara:
	h = 0.704
	OM = 0.27

#DL = distance.luminosity_distance(cosmopara,redshift)*1e6 # pc

table['SB_exp_m'] = -5./2.*(table.SB_Exp+8.9)

## SB cut
#table = table.loc[(table.SB_exp_m < 25),:]

## ---- morphology pick ----- ##

morphology_single = table.loc[(table.Sersic <2) & (table.morphology == 0),:]

table.mor_mix.loc[table.subfindID.isin(morphology_single.subfindID)] = 1 # morphology coverage pick

morphology = table.loc[(table.mor_mix == 1),:]

## ---- kinematics pick ---- ##

disk_star = table.loc[(table.disk_star_f > 0.4),:]
bulge_star = table.loc[(table.bulge_star_f < 0.6),:]
'''
## ---- hybird catagories ---- ##

three = table.loc[(table.mor_mix == 1)&(table.disk_star_f > 0.4)&(table.bulge_star_f < 0.6),:]

#print morphology
#print bulge_star
m_only = morphology.loc[(~morphology.subfindID.isin(bulge_star.subfindID)),:]
k_only = bulge_star.loc[(~bulge_star.subfindID.isin(morphology.subfindID)),:]
b_only = k_only.loc[(~k_only.subfindID.isin(disk_star.subfindID)),:]

three.to_csv(basePath+'redshift_evo/snapshot'+ssNumber+'_three_st.dat',sep = '\t',index = False)
m_only.to_csv(basePath+'redshift_evo/snapshot'+ssNumber+'_mOnly_st.dat',sep = '\t',index = False)
k_only.to_csv(basePath+'redshift_evo/snapshot'+ssNumber+'_kOnly_st.dat',sep = '\t',index = False)
b_only.to_csv(basePath+'redshift_evo/snapshot'+ssNumber+'_bOnly_st.dat',sep = '\t',index = False)

## ---- edge-on & face-on ----- ##
#edge_mor = morphology.loc[(table.theta > 80.) & (table.theta < 100.),:]
#face_mor = morphology.loc[(table.theta < 80.) | (table.theta > 100.),:]
'''
'''
## ---- magnitude cut ---- ##
morphology_cut = morphology.loc[(morphology.r_mag<25.),:]
print morphology.shape, morphology_cut.shape
'''


## ----- numpy hist plot ------ ##

#All = np.array(all_table.R_E)
#Dan = np.array(morphology.R_E)
#df,bf = np.array(disk_star.R_E),np.array(bulge_star.R_E)

All = np.array(all_table.log_stellar_mass)
Dan = np.array(morphology.log_stellar_mass)
df,bf = np.array(disk_star.log_stellar_mass),np.array(bulge_star.log_stellar_mass)

#print All
#print Dan

# histogram
#se = np.linspace(0.0,1.5,15)
se = np.linspace(9,12,20)
#se = np.linspace(50,400,20)
dot = []
for i in range(se.size-1):
	dot.append((se[i]+se[i+1])/2.)

All = np.histogram(All,bins = se)[0].astype(float)
Dan = np.histogram(Dan,bins = se)[0].astype(float)
df = np.histogram(df,bins = se)[0].astype(float)
bf = np.histogram(bf,bins = se)[0].astype(float)

hist = pd.DataFrame({'log_stellar_mass':dot,'All':All,'morphology':Dan,'disk_star':df,'bulge_star':bf})

hist.to_csv(basePath+'redshift_evo/snapshot'+ssNumber+'st_hist_30.dat',sep = '\t',index = False)


'''
plt.plot(dot,Dan/All,color='k',label = 'morphology',linestyle='steps')
plt.plot(dot,bf/All,'b-.',label = 'Bulge Star Fraction < 60%',linestyle = 'steps')
plt.plot(dot,df/All,'r-.',label = 'Disk Star Fraction > 40%',linestyle = 'steps')

plt.xlabel('Einstein Radius (arc sec)')
plt.ylabel('counts')
plt.title('snapshot'+ssNumber+'(z='+str(redshift)+') SB + flux cut')
plt.legend(scatterpoints=1,loc =1)
#plt.legend(loc = 1)
plt.xlim(0.2,1.2)
#plt.ylim(18,32)
#plt.show()
plt.savefig('./illustris_cross/snapshot'+ssNumber+'_hist_dcut.png')
'''

'''
## ----- ploting w/ seaborn ----- ##

fig = sns.jointplot(x = 'V_mag', y = 'color', data=morphology,color = 'b')
fig.set_axis_labels('V','B-V')
#plt.title('snapshot'+ssNumber+'(z='+str(redshift)+')')
#plt.legend(scatterpoints=1,loc =1)
plt.xlim(0,30)
plt.show()
#plt.savefig('./illustris_cross/snapshot'+ssNumber+'_SB.png')
'''

'''
#plt.scatter(bulge_star.log_stellar_mass,bulge_star.fcgs,c = 'b',label = 'Bulge Star Fraction < 60%',marker = '.')
#plt.scatter(morphology.log_stellar_mass,morphology.fcgs,marker = 'x',c='g',label = 'morphology')
#plt.scatter(three.log_stellar_mass,three.fcgs,c = 'r',label = 'satisfy all')
#plt.scatter(face_mor.log_stellar_mass,face_mor.SB_exp_m,c = 'r',marker = 'x',label = 'face-on')
#plt.scatter(edge_mor.log_stellar_mass,edge_mor.SB_exp_m,c = 'b',label = 'edge-on')

plt.scatter(bulge_star.V_mag,bulge_star.color,c = 'b',label = 'Bulge Star Fraction < 60%',marker = '.')
plt.scatter(morphology.V_mag,morphology.color,marker = 'x',c='g',label = 'morphology')
plt.scatter(three.V_mag,three.color,c = 'r',label = 'satisfy all')

plt.title('snapshot'+ssNumber+'(z='+str(redshift)+')')
plt.legend(loc = 1)
#plt.ylim(0,1.2)
plt.xlabel('V_mag')
plt.ylabel('B-V')
#plt.ylim(38,20)
#plt.show()
plt.savefig('./illustris_cross/snapshot'+ssNumber+'_color.png')
'''