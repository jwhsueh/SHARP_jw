import pandas as pd
import seaborn as sns
import numpy as np

import matplotlib.pyplot as plt
import DistanceTool as distance

basePath = '../../data/illustris_1/'
ssNumber = 99
redshift = 0.6

catalog_x = basePath+'full_'+str(ssNumber)+'_x.dat'
catalog_y = basePath+'full_'+str(ssNumber)+'_y.dat'
catalog_z = basePath+'full_'+str(ssNumber)+'_z.dat'
table_x,table_y,table_z = pd.read_csv(catalog_x,sep = '\t'),pd.read_csv(catalog_y,sep = '\t'),pd.read_csv(catalog_z,sep = '\t')

## to a long table
table = table_x.append(table_y,ignore_index = True)
table = table.append(table_z,ignore_index = True)

## relaxation
table = table.loc[(table.relaxation == 0),:]

## velocity dispersion pick
table = table.loc[(table.velDisp > 130) & (table.velDisp < 370),:]

## add log mass
table['log_stellar_mass'] = pd.Series(np.log10(table.stellar_mass))

## add apparent surface brightness
''' cosmopara '''
class cosmopara:
	h = 0.704
	OM = 0.27

DL = distance.luminosity_distance(cosmopara,redshift)*1e6 # pc

table['SB_exp_m'] = pd.Series(table.SB_Exp+5.*(np.log10(DL)-1))

## add disk total flux
disk_M = table.SB_Exp -2.5*np.log10(np.pi*table.SteR_halfmass**2)


table['disk_mag_r'] = pd.Series(disk_M+5.*(np.log10(DL)-1))


## ---- morphology pick ----- ##

morphology = table.loc[(table.Sersic <2) & (table.morphology == 1),:]

## ---- kinematics pick ---- ##

disk_star = table.loc[(table.disk_star_f > 0.4),:]
bulge_star = table.loc[(table.bulge_star_f < 0.6),:]

## ----- ploting w/ seaborn ----- ##

fig = sns.jointplot(x = 'log_stellar_mass', y = 'r_mag', data=morphology,color = 'r')
fig.set_axis_labels('log(stellar_mass)','total r_band mag')
plt.ylim(32,20)
plt.show()