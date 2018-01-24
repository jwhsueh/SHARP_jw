import numpy as np

lens = 'B1422'
fsub = '0533'
path = '/Volumes/sting_1/subs/'+lens+'/result_new2/'+fsub+'/'


idx_list = np.arange(0,50)
idx_list = idx_list.astype(str)

chi2_list = np.array([])

for idx in idx_list:
	filename = path+lens+'_'+fsub+'_'+idx+'out.txt'
	#filename = path+lens+'_abc'+fsub+'_'+idx+'chi2.txt'
	table = np.loadtxt(filename).astype(float)
	#print table.shape
	chi2_list= np.append(chi2_list,table)

chi2_list = np.array(chi2_list)
#print chi2_list.shape
chi2_list = chi2_list.reshape((len(chi2_list)/16,16))
#print chi2_list.shape
np.savetxt(path+lens+'_'+fsub+'_out_com.txt',chi2_list,fmt='%f')