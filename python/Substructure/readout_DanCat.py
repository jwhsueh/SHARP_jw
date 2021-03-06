import numpy as np
import groupcat

ssNumber = '085'
print ssNumber

catPath = '/Volumes/narsil_1/jwhsueh/illustris_1/snapshot_catalog/LensingCataSnap'+str(ssNumber)+'/'
basePath = '../../data/illustris_1/'

Proj = ['NatProj_1/','NatProj_2/','NatProj_3/']

## catalog file names

#galaxy_cata = basePath+'Galaxy_'+str(ssNumber)+'_sig.dat'
#GalaxyID = np.loadtxt(galaxy_cata,dtype = 'int',unpack=True, usecols=[0])

#GroupFirstSub = groupcat.loadHalos(basePath,ssNumber,fields = ['GroupFirstSub'])

TotList = ['TotCata_Lens.list','TotCata_Phot.list','TotCata_Gas.list']
NameList = ['Tot_Lens','Tot_Phot','Tot_Gas']
list_idx = 0
print NameList[list_idx]
proj_name = ['x','y','z']
p_index = 0 # Proj index
print p_index

cataName = basePath+NameList[list_idx]+ssNumber+'_'+proj_name[p_index]+'.dat'

CataList = open(catPath+Proj[p_index]+TotList[list_idx],'r')
CataList = CataList.read().splitlines()

## Lens
field = ['SubfindID','central galaxy = 0','mass [M_sun/h]','R_ein', 'DMfrac w/i R_ein']
#f_idx = np.array([0,1,5,7,11])

## Phot
#field = ['SubfindID','central galaxy = 0','stellar mass','half-stellar mass radius (arcsec)','Sersic index','surface brightness ar Reff_Exp (mag/arcsec^2)','1:Early(deV);0:Late type(Exp)']
#f_idx = np.array([0,1,8,10,16,22,24])

## Gas
#field = ['subfindID','fgas','fcgs']
#f_idx = np.array([0,1,25,28])


new_cata = open(cataName,'w')
#new_cata.write('# '+str(field)+'\n')

for file_name in CataList:

	cata_file = catPath+Proj[p_index]+file_name
	#cata_file.readline() # get rid of the first line

	frs_idx = np.loadtxt(cata_file,dtype = 'int',skiprows = 1,unpack=True, usecols=[1]) # indicate if is the first subhalo in halo

	#cata_file = open(catPath+Proj[p_index]+file_name,'r')
	SubfindID = np.loadtxt(cata_file,dtype = 'int',skiprows=1,unpack=True, usecols=[0])

	#cata_file = open(catPath+Proj[p_index]+file_name,'r')
	sub_table = np.loadtxt(cata_file,skiprows=1,unpack=True, usecols=[1,5,7,11])
	#sub_table = np.loadtxt(cata_file,skiprows=1,unpack=True, usecols=[1,8,10,16,22,24])
	#sub_table = np.loadtxt(cata_file,skiprows=1,unpack=True, usecols=[1,25,28])
	#print sub_table
	valid_id = 2

	if SubfindID.size >1:

		for i in range(len(frs_idx)):
			#if frs_idx[i] == 0 and sub_table[valid_id,i]>0:
			if sub_table[valid_id,i]>0:	

				new_cata.write(str(SubfindID[i]+frs_idx[i])+'	'+str(sub_table[:-2,i])[1:-1]+'	'+str(sub_table[-2:,i])[1:-1]+'\n')


	elif sub_table[valid_id]>0:

		new_cata.write(str(SubfindID+frs_idx)+'	'+str(sub_table[:-2])[1:-1]+' '+str(sub_table[-2:])[1:-1]+'\n')
	
	#cata_file.close()



new_cata.close()

# sort by subfindID
unsort = np.loadtxt(cataName)
sort = unsort[np.argsort(unsort[:,0])]

ids = sort[:,0].astype(int)
print ids[:5]

sort = sort[:,1:]

write_file = open(cataName,'w')

write_file.write('#'+str(field)[1:-1]+'\n')

for i in range(sort.shape[0]):
	write_file.write(str(ids[i])+'	'+str(sort[i,:-2])[1:-1]+' '+str(sort[i,-2:])[1:-1]+'\n')


#for i in 

