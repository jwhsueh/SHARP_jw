import numpy as np
from scipy.stats import multivariate_normal
'''
#B1422

x_data = np.array([0.95065,0.0,0.38925,-0.33388])
y_data = np.array([-0.80215,0.0,0.31998,-0.74771])
err_data = np.array([0.0005,0.0005,0.0005,0.0005])
f_data = np.array([0.024,1.062,1.0,0.551])
ferr_data = np.array([0.006,0.009,1e-5,0.007])
f_idx, p_idx = 2,1
'''
#mock3
x_data = np.array([9.507944e-01,  1.493761e-03,   3.892554e-01, -3.340559e-01])
y_data = np.array([-8.022173e-01 , 1.813538e-03,  3.199820e-01 ,-7.484454e-01])
err_data = np.array([0.0003,0.0003,0.0003,0.0003])
f_data = np.array([0.03412,1.04549,1.0,0.556])
f_data=np.abs(f_data)
f_idx, p_idx = 2,1
f_data = f_data/f_data[f_idx]
ferr_data = f_data*0.03


lens = 'mock3'
fsub = '1000'
path = '/Volumes/sting_1/subs/'+lens+'/result_new3/'

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
print obs_array

cov_table = np.cov(table,rowvar=False)
mid = np.median(table,axis =0)
print mid

## draw for 1000 times
sam = np.zeros(1000)
for i in range(1000):
	obs_set = np.random.multivariate_normal(mid,cov_table,1000)

	#y = multivariate_normal.pdf(obs_array,mean=mid, cov=cov_table)
	obs_cov = np.diag(err_array)
	y = multivariate_normal.pdf(obs_set,mean=mid, cov=obs_cov)
	sam[i] = np.sum(y)/1000

print np.median(sam)
print np.std(sam)

y2= multivariate_normal.pdf(mid,mean=mid, cov=cov_table)
print y2

