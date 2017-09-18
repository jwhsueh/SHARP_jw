from multiprocessing import Pool
import os
from multiprocessing import Process
import numpy as np

def mcmc_run(num):
    """thread function"""
    print 'random number', np.random.rand(1)
    return

#namelist = ['bob','alice','dave','harry','ken']

if __name__ == '__main__':
	jobs = []
	
	for i in range(3):
		#print i
		p = Process(target=mcmc_run, args=(i,))
		#jobs.append(p)
		p.start()
    #p.join()