import numpy as np

def gauss(x,sig):
	return np.exp(-np.power(x, 2.) / (2 * np.power(sig, 2.)))

#lens_list = np.array([0,1,2,3,4,5,6,7,9,10,12,13,14,16,17,18,19,20])
#lens_list = np.array([61,62,63,65,66,67,68,69,70,71,72,74,75,76,77,78,79])
lens_list = np.array([0])

#fsub = '00200sub'
#fsub_list = ['00100_c2','00133_c2','00145_c2','00164_c2','00113_c2']
#fsub_list = ['00100_a','00133_a','00145_a','00164_a','00113_a']
#fsub_list = ['00200_c2','00233_c2','00245_c2','00264_c2','00213_c2']
#fsub_list = ['00200_a','00233_a','00245_a','00264_a','00213_a']
fsub_list = ['00120','00133','00145','00164','00113','00100']
#fsub_list = ['00300_a','00333_a','00345_a','00364_a','00313_a']


#print fsub
#path = '/Volumes/sting_1/subs/mock_test/mc_file/'
path = '/Volumes/narsil_1/jwhsueh/pylens_mc/mock_test/mc_file/'

mock_data = np.loadtxt('/Volumes/narsil_1/jwhsueh/pylens_mc/mock_test/m1_raytrace5.txt')
#print mock_data
mock_x,mock_y,mock_f = mock_data[:,:4],mock_data[:,4:8],mock_data[:,8:12]
smooth_data = np.loadtxt('/Volumes/narsil_1/jwhsueh/pylens_mc/mock_test/m2_raytrace_los.txt')
smooth_x,smooth_y,smooth_f = smooth_data[:,:4],smooth_data[:,4:8],smooth_data[:,8:12]


for lens_idx in lens_list:
	lens = 'mock'+str(lens_idx)
	list_idx = lens_idx
	path1 = path+lens+'/'

	#mock3
	x_data = mock_x[list_idx,:]
	y_data = mock_y[list_idx,:]
	#err_data = np.array([0.003,0.003,0.003,0.003])
	err_data = np.array([0.01,0.01,0.01,0.01])
	f_data = mock_f[list_idx,:]
	f_data=np.abs(f_data)
	f_idx = 1
	f_data = f_data/f_data[f_idx]
	ferr_data = f_data*0.01

	## smooth model

	x_sm,y_sm,f_sm = smooth_x[list_idx,:],smooth_y[list_idx,:],smooth_f[list_idx,:]
	f_sm=np.abs(f_sm)
	f_sm = f_sm/f_sm[f_idx]

	for fsub in fsub_list:
		filename = path+lens+'_'+fsub+'_out_com.txt'
		#filename = path+lens+'_'+fsub+'_0out.txt'
		table = np.genfromtxt(filename)[:,:12]
		table[:,8:] = np.abs(table[:,8:])
		fi_array = np.array([table[:,8+f_idx],table[:,8+f_idx],table[:,8+f_idx],table[:,8+f_idx]]).T
		table[:,8:] = table[:,8:]/fi_array

		#print table[:10,8:]

		obs_array = np.concatenate((x_data,y_data,f_data))
		smooth_array = np.concatenate((x_sm,y_sm,f_sm))
		err_array = np.concatenate((err_data,err_data,ferr_data))

		#print obs_array
		#print smooth_array

		diff_table = table-obs_array
		#print diff_table[:10,8:]

		## position 3-sigma mask

		diff_table = np.abs(diff_table)[:,8:]
		mask_table = -(diff_table-err_array[8:]*3.0)
		#mask_table = -(diff_table-err_array[8:]*10.0)
		#mask_table = -(diff_table-err_array[8:]*2.0)
		#diff_table = np.abs(diff_table)
		#mask_table = -(diff_table-err_array*3.0)

		#print mask_table[:10,:]
		mask_table = mask_table>0

		mask_array = np.all(mask_table,axis=1)
		#print mask_array.shape
		#print np.any(mask_array)

		mask_table = np.repeat([mask_array],12)
		mask_table = mask_table.reshape((len(table[:,0]),12))
		#print mask_table.shape

		table = table[mask_table]
		table = table.reshape((len(table)/12,12))
		#print table.shape

		#print np.min(table[:,9])
		#print table[:10,:]

		diff_table = table-obs_array
		#print diff_table[:10,:8]

		## position gauss prior

		pos_pri = np.exp(-0.5*(diff_table[:,:8])**2/(err_array[:8])**2)
		pos_pri = np.prod(pos_pri,axis=1)
		pos_flat = np.exp(-0.5*3.0**2)

		#pos_pri = np.exp(-0.5*(diff_table)**2/(err_array)**2)
		#pos_pri = np.prod(pos_pri,axis=1)
		#print pos_pri[:10]

		diff_table = diff_table**2
		err_array=err_array**2

		#chi2_table = diff_table[:,:8]/err_array[:8]
		chi2_table = diff_table[:,8:]/err_array[8:]
		#chi2_table = diff_table/err_array
		chi2_list = np.sum(chi2_table,axis=1)/11.
		#print chi2_list.shape
		
		#idx = np.argmin(chi2_list)
		#chi2_list=np.delete(chi2_list,idx)
		#pos_pri=np.delete(pos_pri,idx)
		'''
		### smooth model likelihood
		smooth_chi2 = (smooth_array-obs_array)**2/err_array
		smooth_chi2 = smooth_chi2[8:]
		smooth_chi2 = np.sum(smooth_chi2)/11.
		#print smooth_chi2
		smooth_likehood = np.exp(-smooth_chi2/2.0)
		#print smooth_likehood

		likehood = np.exp(-chi2_list/2.0)
		#likehood = likehood-smooth_likehood

		mask2 = likehood>smooth_likehood
		likehood = likehood[mask2]
		pos_pri = pos_pri[mask2]
		'''
		#likehood = np.exp(-chi2_list/2.0)*pos_pri
		#likehood[mask_array] = np.exp(-chi2_list[mask_array]/2.0)*pos_flat
		#likehood = np.exp(-chi2_list[mask_array]/2.0)

		likehood = np.exp(-chi2_list/2.0)*pos_pri
		#likehood = likehood#*pos_pri

		#likehood = pos_pri

		#np.savetxt(path+lens+'_'+fsub+'_pri_likelihood.txt',np.c_[chi2_list,likehood])
		#np.savetxt(path+lens+'_'+fsub+'_pri_table.txt',table)
		#print np.min(chi2_list)
		#print len(likehood)
		#print likehood[:10]
		#print np.sum(likehood)
		print np.sum(likehood)*len(likehood)
	#print np.sum(likehood)/len(likehood)
	#print np.sum(likehood)/len(likehood)*np.exp(-np.min(chi2_list)/2.0)
