import numpy as np

path = '/Volumes/sting_1/subs/B1422/protect/'

idx_list = np.arange(3,15)
idx_list = idx_list.astype(str)

chi2_list = np.array([])

for idx in idx_list:
	filename = path+'B1422_abc0200_'+idx+'chi2.txt'
	table = np.loadtxt(filename).astype(float)
	chi2_list= np.append(chi2_list,table)

chi2_list = np.array(chi2_list)
np.savetxt(path+'B1422_abc0200_chi2_com2.txt',chi2_list,fmt='%f')