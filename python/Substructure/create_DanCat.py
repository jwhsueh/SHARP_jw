import pandas as pd
import numpy as np

## read-in Dandan's catalog list and create two 'TotList'

ssNumber = '120'

catPath = '/Volumes/narsil_1/jwhsueh/illustris_1/snapshot_catalog/LensingCataSnap'+ssNumber+'/'

catfile = 'InAll3Proj.list'

table = pd.read_csv(catPath+catfile,sep= '_',names = 'abcdef')

# get the list of fxx
# get rid of duplicates
table = table.drop_duplicates('c')

#print table

name_list = np.array(table.c)

## output Lens & Phot file list

TotList = ['TotCata_Lens.list','TotCata_Phot.list','TotCata_Gas.list']

LensList = open(catPath+TotList[0],'w')
PhotList = open(catPath+TotList[1],'w')
GasList = open(catPath+TotList[2],'w')

for f_name in name_list:
	LensList.write('IllustrisLens_'+f_name+'_LensingGalParam.dat\n')
	PhotList.write('IllustrisLens_'+f_name+'_LightMeasInStrDymBID.dat\n')
	GasList.write('IllustrisLens_'+f_name+'_ShapeFracSlpMeasInMorphBID.dat\n')

LensList.close()
PhotList.close()
GasList.close()