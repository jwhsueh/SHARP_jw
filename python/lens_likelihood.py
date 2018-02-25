import numpy as np
'''
# B1422

x_data = np.array([0.95065,0.0,0.38925,-0.33388])
y_data = np.array([-0.80215,0.0,0.31998,-0.74771])
#err_data = np.array([0.0005,0.0005,0.0005,0.0005])
err_data = np.array([0.005,0.005,0.005,0.005])*0.6
f_data = np.array([0.024,1.062,1.0,0.551])
ferr_data = np.array([0.006,0.009,1e-5,0.007])
#ferr_data = f_data*0.03
f_idx, p_idx = 2,1
'''

lens_idx = 0
lens = 'B1422'
fsub = '1000'
path = '/Volumes/sting_1/subs/lens_info/'
path1 = '/Volumes/sting_1/subs/'+lens+'/result_new4/'

mock_data = np.loadtxt(path+'lens_obs.txt')
mock_x,mock_y,mock_f,mock_ferr = mock_data[:,:4],mock_data[:,4:8],mock_data[:,8:12],mock_data[:,12:16]

lens_fidx = np.loadtxt(path+'lens_fidx.txt')

#mock3
x_data = mock_x[lens_idx,:]
y_data = mock_y[lens_idx,:]
#err_data = np.array([0.003,0.003,0.003,0.003])
err_data = np.array([0.01,0.01,0.01,0.01])
f_data = mock_f[lens_idx,:]
f_data=np.abs(f_data)
f_idx = lens_fidx[lens_idx]
f_data = f_data/f_data[f_idx]
ferr_data = mock_ferr[lens_idx,:]



filename = path1+lens+'_'+fsub+'_out_com.txt'
#filename = path+lens+'_'+fsub+'_0out.txt'
table = np.genfromtxt(filename)[:,:12]
table[:,8:] = np.abs(table[:,8:])
fi_array = np.array([table[:,8+f_idx],table[:,8+f_idx],table[:,8+f_idx],table[:,8+f_idx]]).T
table[:,8:] = table[:,8:]/fi_array


#print table[:10,8:]

obs_array = np.concatenate((x_data,y_data,f_data))
err_array = np.concatenate((err_data,err_data,ferr_data))

diff_table = table-obs_array
#print diff_table[:10,8:]

## position 3-sigma mask

def gauss(x,sig):
	return np.exp(-np.power(x, 2.) / (2 * np.power(sig, 2.)))

diff_table = np.abs(diff_table)[:,8:]
mask_table = -(diff_table-err_array[8:]*3.0)

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
'''
idx = np.argmin(chi2_list)
chi2_list=np.delete(chi2_list,idx)
pos_pri=np.delete(pos_pri,idx)
'''

#likehood = np.exp(-chi2_list/2.0)*pos_pri
#likehood[mask_array] = np.exp(-chi2_list[mask_array]/2.0)*pos_flat
#likehood = np.exp(-chi2_list[mask_array]/2.0)
likehood = np.exp(-chi2_list/2.0)*pos_pri
#likehood = pos_pri

#np.savetxt(path+lens+'_'+fsub+'_pri_likelihood.txt',np.c_[chi2_list,likehood])
#np.savetxt(path+lens+'_'+fsub+'_pri_table.txt',table)
print np.min(chi2_list)
print len(likehood)
#print likehood[:10]
print np.sum(likehood)*len(likehood)