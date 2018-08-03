import numpy as np
import os

lens = 'mock9'
#lens = 'MG0414'
fsub = '0200s'
path = './'
path1 = path
#path = '/Volumes/narsil_1/jwhsueh/pylens_mc/mock_test/mc_file/'
#path1 = path+lens+'/'


idx_list = np.arange(0,40)
idx_list = idx_list.astype(str)

chi2_list = np.array([])
step_list = 0

#### ----- get rid of last line ------ ####
for idx in idx_list:
	txt = open(path1+lens+'_'+fsub+'_'+idx+'out.txt','r')
	lines = txt.readlines()
	lines = lines[:-1]
	txt.close()

	txt = open(path1+lens+'_'+fsub+'_'+idx+'out_2.txt','w')
	txt.writelines(lines)
	txt.close()

for idx in idx_list:
	print idx
	filename = path1+lens+'_'+fsub+'_'+idx+'out_2.txt'
	#filename = path1+lens+'_abc'+fsub+'_'+idx+'chi2.txt'
	#filename2 = path1+lens+'_'+fsub+'_'+idx+'step.txt'
	table = np.genfromtxt(filename).astype(float)
	#step = np.loadtxt(filename2)
	#print table.shape
	chi2_list= np.append(chi2_list,table)
	#step_list = step_list+step

chi2_list = np.array(chi2_list)
#print step_list
chi2_list = chi2_list.reshape((len(chi2_list)/16,16))
print chi2_list.shape[0]
np.savetxt(path+lens+'_'+fsub+'_out_com.txt',chi2_list,fmt='%f')
os.system('rm '+path1+lens+'_'+fsub+'_'+'*out_2.txt ')
os.system('rm '+path1+lens+'_'+fsub+'_'+'*step.txt ')