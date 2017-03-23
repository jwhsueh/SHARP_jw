import numpy as np
import matplotlib.pyplot as plt

n_src=50

path='/Volumes/sting_1/data/'
datapath=path+'cm_test/'

list_file = path+'snap99_tri2_Rcusp.txt'

file_list = np.genfromtxt(list_file,dtype='str')

for file_one in file_list:
	file_str = file_one.split('_p')
	subfindID = file_str[0]
	print subfindID

	caus=np.loadtxt('/Volumes/sting_1/snap99_'+subfindID+'/caustic_'+file_one+'.txt')
	#caus2=np.loadtxt('/Volumes/sting_1/snap99_275833/caustic_275833_p1sub.txt')
	#crit=np.loadtxt('/Users/jwhsueh/Documents/glamer/examples/ParticleExample/build/critical_proj3_281185sub_64.txt')
	#crit2=np.loadtxt('/Users/jwhsueh/Documents/glamer/examples/ParticleExample/build/critical_proj3_281185_64.txt')
	src_pt=np.loadtxt('/Volumes/sting_1/snap99_'+subfindID+'/'+file_one+'_src.dat')
	#src2_pt=np.loadtxt('/Volumes/sting_1/snap99_275833/275833_p1sub_src.dat')
	drop = np.loadtxt('/Volumes/sting_1/snap99_'+subfindID+'/'+file_one+'_64_drop.txt')
	drop = np.intersect1d(drop,drop)

	#mask2=crit2[:,2]==0
	#mask=crit[:,2]==0
	#caus_intv=np.sqrt((caus[1:,0]-caus[:-1,0])**2+(caus[1:,1]-caus[:-1,1])**2)
	#print caus[1:,0]-caus[:-1,0]
	#print min(caus_intv)

	## ----- calculating distance to caustic

	dist_caus = np.zeros(len(src_pt[:,0]))
	dist = np.zeros(len(caus[:,0]))

	for i in range(len(dist_caus)):
		sx,sy = src_pt[i,0],src_pt[i,1]
		for j in range(len(caus[:,0])):
			cax,cay = caus[j,0],caus[j,1]
			dist[j] = np.sqrt((sx-cax)**2+(sy-cay)**2)

		dist_caus[i] = np.min(dist[j])

	print np.max(dist_caus),np.min(dist_caus)

	hf_dist = (np.max(dist_caus)+np.min(dist_caus))/2
	mid_dist = np.median(dist_caus)

	hf_mask = dist_caus<hf_dist
	mid_mask = dist_caus<mid_dist

	full_mask = np.arange(n_src)
	drop_mask = ~(np.in1d(full_mask,drop))

	hf_mask = hf_mask[drop_mask]
	mid_mask = mid_mask[drop_mask]

	print hf_mask.size, mid_mask.size

	np.savetxt(datapath+'snap99_'+file_one+'_hfmask.txt',hf_mask,fmt='%d')
	np.savetxt(datapath+'snap99_'+file_one+'_midmask.txt',mid_mask,fmt='%d')

