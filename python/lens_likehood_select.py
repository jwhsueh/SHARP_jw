import numpy as np

lens_idx = 3
lens = 'PG1115'
fsub = '0200'
path = '/Volumes/narsil_1/jwhsueh/pylens_mc/lens_info/'
path1 = '/Volumes/narsil_1/jwhsueh/pylens_mc/'+lens+'/'

select_N = 90000

mock_data = np.loadtxt(path+'lens_obs.txt')
mock_x,mock_y,mock_f,mock_ferr = mock_data[:,:4],mock_data[:,4:8],mock_data[:,8:12],mock_data[:,12:16]

lens_fidx = np.loadtxt(path+'lens_fidx.txt')

#mock3
x_data = mock_x[lens_idx,:]
y_data = mock_y[lens_idx,:]
#err_data = np.array([0.003,0.003,0.003,0.003])
err_data = np.array([0.01,0.01,0.01,0.01])
#err_data = np.array([0.1,0.1,0.1,0.1])
f_data = mock_f[lens_idx,:]
f_data=np.abs(f_data)
f_idx = lens_fidx[lens_idx]
f_data = f_data/f_data[f_idx]
ferr_data = mock_ferr[lens_idx,:]
#print f_data

filename = path1+lens+'_'+fsub+'_out_com.txt'
#filename = path+lens+'_'+fsub+'_0out.txt'
table = np.genfromtxt(filename)[:,:12]
table[:,8:] = np.abs(table[:,8:])
fi_array = np.array([table[:,8+f_idx],table[:,8+f_idx],table[:,8+f_idx],table[:,8+f_idx]]).T
table[:,8:] = table[:,8:]/fi_array
idx = np.arange(len(table[:,0]))

#print table[:10,8:]

obs_array = np.concatenate((x_data,y_data,f_data))
err_array = np.concatenate((err_data,err_data,ferr_data))

def gauss(x,sig):
	return np.exp(-np.power(x, 2.) / (2 * np.power(sig, 2.)))


#####

## position 3-sigma mask
diff_table = table-obs_array
diff_table = np.abs(diff_table)[:,8:]
mask_table = -(diff_table-err_array[8:]*3.0)
mask_table = mask_table>0
mask_array = np.all(mask_table,axis=1)
print mask_array.shape

mask_table = np.repeat([mask_array],12)
mask_table = mask_table.reshape((len(table[:,0]),12))

table = table[mask_table]
table = table.reshape((len(table)/12,12))

diff_table = table-obs_array

## position gauss prior

pos_pri = np.exp(-0.5*(diff_table[:,:8])**2/(err_array[:8])**2)
pos_pri = np.prod(pos_pri,axis=1)
pos_flat = np.exp(-0.5*3.0**2)

diff_table = diff_table**2
err_array=err_array**2

#chi2_table = diff_table[:,:8]/err_array[:8]
chi2_table = diff_table[:,8:]/err_array[8:]
#chi2_table = diff_table/err_array
chi2_list = np.sum(chi2_table,axis=1)/11.
likehood = np.exp(-chi2_list/2.0)*pos_pri

#### selection
N = 100
likehood_array=np.zeros(N)

chi2_id = idx[mask_array]

for i in range(N):

	select_id = np.random.choice(idx,select_N)
	select_mask = np.in1d(chi2_id,select_id)

	likehood_select = likehood[select_mask]
	likehood_array[i] = np.sum(likehood_select)*len(likehood_select)

print np.sort(likehood_array)
print np.average(likehood_array)
print np.median(likehood_array)