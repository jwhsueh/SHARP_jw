import numpy as np

#B1422

x_data = np.array([0.95065,0.0,0.38925,-0.33388])
y_data = np.array([-0.80215,0.0,0.31998,-0.74771])
err_data = np.array([0.0005,0.0005,0.0005,0.0005])
f_data = np.array([0.024,1.062,1.0,0.551])
ferr_data = np.array([0.006,0.009,1e-5,0.007])
f_idx, p_idx = 2,1

lens = 'B1422'
fsub = '00200'
path = '/Volumes/sting_1/subs/'+lens+'/result_new2/'

filename = path+lens+'_'+fsub+'_out_com.txt'
table = np.loadtxt(filename)[:,:12]
table[:,8:] = np.abs(table[:,8:])
fi_array = np.array([table[:,8+f_idx],table[:,8+f_idx],table[:,8+f_idx],table[:,8+f_idx]]).T
table[:,8:] = table[:,8:]/fi_array

obs_array = np.concatenate((x_data,y_data,f_data))
err_array = np.concatenate((err_data,err_data,ferr_data))

table = np.delete(table,8+f_idx,axis=1)
obs_array = np.delete(obs_array,8+f_idx)
err_array = np.delete(err_array,8+f_idx)

diff_table = np.abs(table-obs_array)/err_array
diff_table = np.sum(diff_table,axis=0)
diff_table = 1./diff_table

print np.sum(diff_table)

