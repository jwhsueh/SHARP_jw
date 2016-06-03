import numpy as np

catPath = '../../data/Snap099Data4JW/Snap099Cata/'
basePath = '../../data/illustris_1/'

ssNumber = 99

Proj = ['NatProj_1/','NatProj_2/','NatProj_3/']

## catalog file names

galaxy_cata = basePath+'Galaxy_'+str(ssNumber)+'.dat'
GalaxyID = np.loadtxt(galaxy_cata,dtype = 'int',unpack=True, usecols=[0])

TotList = ['TotCata_Lens.list','TotCata_Phot.list']

p_index = 0 # Proj index

CataList = open(catPath+Proj[p_index]+TotList[1],'r')
CataList = CataList.read().splitlines()

field = ['SubfindID','mass','R_ein', 'DMfrac w/i R_ein']
f_idx = np.array([0,6,7,11])

new_cata = open(basePath+'Dandan_Phot'+str(ssNumber)+'.dat','w')
#new_cata.write('# '+str(field)+'\n')

for file_name in CataList:

	cata_file = open(catPath+Proj[p_index]+file_name,'r')
	#cata_file.readline() # get rid of the first line

	cataID = np.loadtxt(cata_file,dtype = 'int',skiprows = 1,unpack=True, usecols=[0])
	print cataID

	cata_file = open(catPath+Proj[p_index]+file_name,'r')
	sub_table = np.loadtxt(cata_file,skiprows=1,unpack=True, usecols=[0,6,7,11])
	#print cataID.size
	#print sub_table.shape
	
	if cataID.size>1:
		for i in range(len(cataID)):
			if cataID[i] in GalaxyID:
				#print cataID[i],i
				new_cata.write(str(sub_table[:,i])+'\n')

	else:
		if cataID in GalaxyID:
			new_cata.write(str(sub_table)+'\n')
	
	cata_file.close()



new_cata.close()

#for i in 

