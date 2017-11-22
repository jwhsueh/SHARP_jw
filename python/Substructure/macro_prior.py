import numpy as np
from scipy import stats

## lens macro model prior

lens = 'MG0414'

## ---- importance sampling prior

# B1422
#paras = np.array([7.466390e-01, 7.645506e-01, -6.730964e-01, 3.782699e-01, 5.556224e+01, 1.473014e-01, 5.249086e+01,3.902577e-01, -4.156271e-01])
#sig = np.array([0.04,0.02,0.02,0.06,0.73,0.02,1.07,0.02,0.02])

#MG0414
paras = np.array([1.113297e+00, -5.569528e-01, -1.349765e+00, 3.967119e-01, -7.256150e+01, 3.455577e-02, 6.986017e+01,-3.507265e-01, -1.179662e+00])
#sig = np.array([0.04,0.02,0.02,0.06,0.73,0.02,1.07,0.02,0.02])

## draw macro model parameters

N = 100000

# generate covariance matrix from MCMC chain
table = np.loadtxt('../../data/sub_gravlens/'+lens+'_mcmc_param_lens.txt')

cova = np.cov(table,rowvar=False)
mean = np.mean(table,axis =0)

#min_eig = np.min(np.real(np.linalg.eigvals(cova)))
#if min_eig < 0:
#    cova -= 10*min_eig * np.eye(*cova.shape)

para_table = np.zeros((N,len(paras)+1))

for i in range(N):
	para_table[i,:-1] = np.random.multivariate_normal(mean,cova)

for i in range(N):
	macro_para = para_table[i,:-1]
	weight = stats.multivariate_normal.pdf(macro_para,mean,cov=1) ##
	para_table[i,-1] = weight


print para_table[:,-1]

np.savetxt('/Volumes/sting_1/subs/'+lens+'_eta.txt',para_table)