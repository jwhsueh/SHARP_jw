import numpy as np
import groupcat

catPath = '../../data/Snap099Data4JW/Snap099Cata/'
basePath = '../../data/illustris_1/'

ssNumber = 99

Proj = ['NatProj_1/','NatProj_2/','NatProj_3/']

## catalog file names

galaxy_cata = basePath+'Galaxy_'+str(ssNumber)+'.dat'
GalaxyID = np.loadtxt(galaxy_cata,dtype = 'int',unpack=True, usecols=[0])

#GroupFirstSub = groupcat.loadHalos(basePath,ssNumber,fields = ['GroupFirstSub'])

TotList = ['TotCata_Lens.list','TotCata_Phot.list']
list_idx = 0

p_index = 0 # Proj index

CataList = open(catPath+Proj[p_index]+TotList[list_idx],'r')
CataList = CataList.read().splitlines()

## Lens
field = ['SubfindID','mass [M_sun/h]','mass w/i R_ein','R_ein', 'DMfrac w/i R_ein']
f_idx = np.array([0,5,6,7,11])

## Phot
#field = ['SubfindID','mass [M_sun/h]','stellar mass','Sersic index', '1:Early(deV);0:Late type(Exp)']
#f_idx = np.array([0,3,4,18,24])


new_cata = open(basePath+'Dandan_Lens'+str(ssNumber)+'.dat','w')
#new_cata.write('# '+str(field)+'\n')

for file_name in CataList:

	cata_file = catPath+Proj[p_index]+file_name
	#cata_file.readline() # get rid of the first line

	frs_idx = np.loadtxt(cata_file,dtype = 'int',skiprows = 1,unpack=True, usecols=[1]) # indicate if is the first subhalo in halo

	#cata_file = open(catPath+Proj[p_index]+file_name,'r')
	groupID = np.loadtxt(cata_file,dtype = 'int',skiprows = 1,unpack=True, usecols=[0])
	print groupID

	#cata_file = open(catPath+Proj[p_index]+file_name,'r')
	SubfindID = np.loadtxt(cata_file,dtype = 'int',skiprows=1,unpack=True, usecols=[0])

	#cata_file = open(catPath+Proj[p_index]+file_name,'r')
	sub_table = np.loadtxt(cata_file,skiprows=1,unpack=True, usecols=[5,6,7,11])
	valid_id = 1

	if groupID.size >1:

		for i in range(len(frs_idx)):
			if frs_idx[i] == 0 and sub_table[valid_id,i]>0:			

				new_cata.write(str(SubfindID[i])+'	'+str(sub_table[:,i])[1:-1]+'\n')

	elif sub_table[valid_id]>0:
		
		new_cata.write(str(SubfindID)+'	'+str(sub_table)[1:-1]+'\n')
	
	#cata_file.close()



new_cata.close()

sort_file = basePath+'Dandan_Lens'+str(ssNumber)+'.dat'
unsort = np.loadtxt(sort_file)
sort = unsort[np.argsort(unsort[:,0])]

ids = sort[:,0].astype(int)
print ids[:5]

sort = sort[:,1:]

write_file = open(basePath+'Dandan_Lens'+str(ssNumber)+'.dat','w')

write_file.write('#'+str(field)[1:-1]+'\n')

for i in range(sort.shape[0]):
	write_file.write(str(ids[i])+'	'+str(sort[i,:])[1:-1]+'\n')


#for i in 
