import numpy as np
from scipy import stats

## lens macro model prior

lens = 'B1422'

## ---- importance sampling prior
paras = np.array([7.466390e-01, 7.645506e-01, -6.730964e-01, 3.782699e-01, 5.556224e+01, 1.473014e-01, 5.249086e+01,3.902577e-01, -4.156271e-01])
sig = np.array([0.04,0.02,0.02,0.06,0.73,0.02,1.07,0.02,0.02])

## ---- for starting point
paras2 = np.array([7.466390e-01, 7.645506e-01, -6.730964e-01, 3.782699e-01, 5.556224e+01, 1.473014e-01, 5.249086e+01,3.902577e-01, -4.156271e-01])
sig2 = np.array([0.0001,0.0001,0.001,0.001,0.1,0.0001,0.1,0.0001,0.0001])



## draw macro model parameters

N = 110000

# generate covariance matrix from MCMC chain
table = np.loadtxt('../../data/sub_gravlens/B1422_flatchain_0.txt')

cova = np.cov(table,rowvar=False)
mean = np.mean(table,axis =0)

para_table = np.zeros((N,len(paras)+1))

for i in range(N):
	para_table[i,:-1] = np.random.multivariate_normal(mean,cova)

for i in range(N):
	macro_para = para_table[i,:-1]
	weight = stats.multivariate_normal.pdf(macro_para,mean,cov=1) ##
	para_table[i,-1] = weight


print para_table[:,-1]

np.savetxt('/Volumes/sting_1/subs/B1422_eta.txt',para_table)