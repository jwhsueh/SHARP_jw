import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

path = '/Volumes/narsil_1/jwhsueh/pylens_mc/'
lens_list = np.genfromtxt(path+'lens_info/lens_list.txt',dtype=str)
#lens_list = np.array(lens_list[0:2])
#print lens_list

def gaus(x,a,x0,sigma):

	return a*np.exp(-(x-x0)**2/(2*sigma**2))

likelihood = np.zeros(100)
likelihood.fill(1.0)

## colors
colors = ['b','r','m','k','g']

for i in range(len(lens_list)):
	lens = lens_list[i]
	pdf_file = path+lens+'/'+lens+'_pdf.txt'
	like_table = np.loadtxt(pdf_file)

	f_sub_l = like_table[:,0]
	f_sub = np.log10(f_sub_l)
	#f_dot = np.logspace(np.log10(0.005),np.log10(0.02),100)
	f_dot = np.linspace(np.log10(0.002),np.log10(0.1),100)
	f_dot_l = 10**f_dot

	n = float(len(f_sub))                      
	#drop = np.array([])
	#print popt

	#plt.semilogx(f_sub_l,case,marker='^',linestyle='none')
	#plt.semilogx(f_dot_l,gaus(f_dot,*popt),'r:')
	#plt.legend()
	#plt.show()
	print lens
	case = like_table[:,2]
	case = case/np.sum(case)
	mean = np.sum(f_sub*case)/n                   
	sigma = np.sum(case*(f_sub-mean)**2)/n        
	p0=[1.0,mean,sigma]
	plt.semilogx(f_sub_l,case,marker='^',linestyle='-',label=lens,color=colors[i])
	popt,pcov = curve_fit(gaus,f_sub,case,p0)
	plt.semilogx(f_dot_l,gaus(f_dot,*popt),':',color=colors[i])
	likelihood = likelihood*gaus(f_dot,*popt)



'''
joint = np.zeros(3)
for i in range(like_table.shape[1]):
	joint[i] = np.prod(likelihood[:,i])

print joint
joint = joint/np.sum(joint)
'''

## likelihood normalize
print likelihood[50:60]
likelihood = likelihood/np.sum(likelihood)/0.01*(0.1-0.002)

#likelihood2 = likelihood2/np.sum(likelihood2)

#plt.semilogx(f_sub,joint)
plt.xlim(0.001,0.2)
#plt.semilogx(f_dot_l,likelihood)
#plt.semilogx(f_sub_l,likelihood2)
#plt.plot([0.01,0.01],[0,0.3])
plt.xlabel('fsub')
plt.ylabel('normalized likelihood')
plt.title('pdf of real lenses')
plt.legend(scatterpoints=1)
plt.show()
#plt.savefig(path+'lens_info/lens_pdf_fit.png')
