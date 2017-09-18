import numpy as np
import os

N = 5000
sec = 10 # number of jobs
n_sec = int(N/sec)
a = np.linspace(0,N,sec+1).astype(int)
lens='B1422'

for i in range(len(a)-1):
	filename = './submit'
	submit_file = open(filename,'w')
	submit_file.write('#!/bin/bash -l \n\n')
	submit_file.write('#SBATCH --job-name='+lens+'_'+str(i)+'\n')
	submit_file.write('#SBATCH --partition=med\n')
	submit_file.write('#SBATCH --time=24:00:00\n')
	submit_file.write('#SBATCH --nodes=1\n')
	submit_file.write('#SBATCH --mem-per-cpu=2000\n')
	submit_file.write('#SBATCH --ntasks=1\n\n')
	submit_file.write('python -u create_realNFW_mf.py 0.01 01 '+str(n_sec)+' '+str(a[i])+' > '+lens+'_'+str(i)+'.log\n')

	submit_file.close()
	os.system('sbatch submit')