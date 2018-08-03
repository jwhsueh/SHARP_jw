import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

path = '/Volumes/narsil_1/jwhsueh/pylens_mc/mock_test/'
like_table = np.loadtxt(path+'mc_likelihood_set2_sub.txt')

f_sub_l = np.array([0.002,0.005,0.01,0.02,0.05])
f_sub = np.log10(f_sub_l)
#f_dot = np.logspace(np.log10(0.005),np.log10(0.02),100)
f_dot = np.linspace(np.log10(0.002),np.log10(0.05),100)
f_dot_l = 10**f_dot


def gaus(x,a,x0,sigma):

	return a*np.exp(-(x-x0)**2/(2*sigma**2))

n = float(len(f_sub))                      
drop = np.array([8,12,13])
#drop = np.array([])
#print popt

#plt.semilogx(f_sub_l,case,marker='^',linestyle='none')
#plt.semilogx(f_dot_l,gaus(f_dot,*popt),'r:')
#plt.legend()
#plt.show()
likelihood = np.zeros(100)
likelihood.fill(1.0)
likelihood2 = np.zeros(5)
likelihood2.fill(1.0)
for i in range(like_table.shape[0]):
	print i

	
	case = like_table[i,:]
	case = case/np.sum(case)
	mean = np.sum(f_sub*case)/n                   
	sigma = np.sum(case*(f_sub-mean)**2)/n        
	p0=[1.0,mean,sigma]
	#likelihood2 = likelihood2*case
	#plt.semilogx(f_sub_l,case,marker='^',linestyle='-')
	
	if ~np.in1d(i,drop):
		
		popt,pcov = curve_fit(gaus,f_sub,case,p0)
		
		#plt.semilogx(f_dot_l,gaus(f_dot,*popt),':')
		likelihood = likelihood*gaus(f_dot,*popt)
		#likelihood2 = likelihood2*case
	#else:
	#	for j in range(len(f_idx)):
	#		likelihood[f_idx[j]] = likelihood[f_idx[j]]*case[j]

'''
joint = np.zeros(5)
for i in range(like_table.shape[1]):
	joint[i] = np.prod(likelihood[:,i])

print joint
joint = joint/np.sum(joint)
'''

## likelihood normalize
#print likelihood[50:60]
likelihood = likelihood/np.sum(likelihood)/0.01*(0.05-0.002)

## fit to likelihood
mean = np.sum(f_dot*likelihood)/n                   
sigma = np.sum(likelihood*(f_dot-mean)**2)/n        
p0=[1.0,mean,sigma]
popt,pcov = curve_fit(gaus,f_dot,likelihood,p0)
print popt
print 10**popt[1], 10**(popt[1]+popt[2])

#likelihood2 = likelihood2/np.sum(likelihood2)

#plt.semilogx(f_sub,joint)
#plt.xlim(0.002,0.05)
plt.semilogx(f_dot_l,likelihood)
#plt.semilogx(f_sub_l,likelihood2,marker='o')
plt.plot([0.01,0.01],[0,0.3])
plt.ylim(0,0.3)
plt.xlabel('fsub')
plt.ylabel('normalized likelihood')
plt.title('pdf of 20 mock systems [set 2] (fsub=1%, flux err=1%)')
#plt.show()
plt.savefig(path+'likelihood_curve_sub_set2.png')
