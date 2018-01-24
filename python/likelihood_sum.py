import numpy as np

# B1422

x_data = np.array([0.95065,0.0,0.38925,-0.33388])
y_data = np.array([-0.80215,0.0,0.31998,-0.74771])
err_data = np.array([0.0005,0.0005,0.0005,0.0005])
f_data = np.array([0.024,1.062,1.0,0.551])
ferr_data = np.array([0.006,0.009,1e-5,0.007])
f_idx, p_idx = 2,1


lens = 'B1422'
fsub = '0500'
path = '/Volumes/sting_1/subs/'+lens+'/result_new/'

filename = path+lens+'_'+fsub+'_out_com.txt'
#filename = path+lens+'_'+fsub+'_0out.txt'
table = np.loadtxt(filename)[:,:12]
table[:,8:] = np.abs(table[:,8:])
fi_array = np.array([table[:,8+f_idx],table[:,8+f_idx],table[:,8+f_idx],table[:,8+f_idx]]).T
table[:,8:] = table[:,8:]/fi_array

#print table[:10,8:]

obs_array = np.concatenate((x_data,y_data,f_data))
err_array = np.concatenate((err_data,err_data,ferr_data))
err_array=err_array**2

diff_table = table-obs_array
print diff_table[:10,8:]
diff_table = diff_table**2

chi2_table = diff_table/err_array
chi2_list = np.sum(chi2_table,axis=1)/11.
print chi2_list.shape

likehood = np.exp(-chi2_list/2.0)

np.savetxt(path+lens+'_'+fsub+'_likelihood.txt',np.c_[chi2_list,likehood])
print np.min(chi2_list)
print len(likehood)
print likehood[:10]
print np.sum(likehood)