import numpy as np
import matplotlib.pyplot as plt
#from scipy.optimize import curve_fit

path = '/Volumes/narsil_1/jwhsueh/pylens_mc/'
lens_list = np.genfromtxt(path+'lens_info/lens_list.txt',dtype=str)
lens_list = np.array(['B1422','B0128','PG1115','MG0414','Q2237'])
#lens = 'B1422'
#print lens_list

#def gaus(x,a,x0,sigma):
#	return a*np.exp(-(x-x0)**2/(2*sigma**2))

#likelihood = np.zeros(100)
#likelihood.fill(1.0)

## colors
colors = ['b','r','m','k','g']
mass = np.array([1.165e6,1.54e7,5e7,1.4e8,7.4e8])
#mass = np.array([1.4e8,5e7,1.54e7,0.0])
#wdm = np.array([3.3,4.5,6.4,13.9])
mass = np.log10(mass)
f_sub = np.array([0.002,0.005,0.01])
#f_sub = np.array([0.005,0.01,0.015,0.02,0.03,0.05,0.1])
f_sub = np.log10(f_sub)
#wdm = np.array([6.4,13.9])

pdf_com = np.zeros((3,5))
pdf_com.fill(1.0)

for lens in lens_list:
	pdf_file = path+lens+'/'+lens+'_wdm.txt'
	like_table = np.loadtxt(pdf_file)#[:,:-1]
	pdf = like_table/np.sum(like_table)
	pdf_com = pdf_com*pdf


pdf_com = pdf_com/np.sum(pdf_com)
'''
for mass in wdm:
	mask = like_table[:,1] == mass
	f_sub_l = like_table[mask,0]
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
	#print lens
	case = like_table[mask,2]
	case = case/np.sum(case)
	mean = np.sum(f_sub*case)/n                   
	sigma = np.sum(case*(f_sub-mean)**2)/n        
	p0=[1.0,mean,sigma]
	plt.semilogx(f_sub_l,case,marker='^',linestyle='-',label=mass)
	popt,pcov = curve_fit(gaus,f_sub,case,p0)
	plt.semilogx(f_dot_l,gaus(f_dot,*popt),':')
	likelihood = likelihood*gaus(f_dot,*popt)
'''


'''
joint = np.zeros(3)
for i in range(like_table.shape[1]):
	joint[i] = np.prod(likelihood[:,i])

print joint
joint = joint/np.sum(joint)
'''

## likelihood normalize
#print likelihood[50:60]
#likelihood = likelihood/np.sum(likelihood)/0.01*(0.1-0.002)

#likelihood2 = likelihood2/np.sum(likelihood2)

#plt.semilogx(f_sub,joint)
#plt.xlim(0.001,0.2)
#plt.semilogx(f_dot_l,likelihood)
#plt.semilogx(f_sub_l,likelihood2)
#plt.plot([0.01,0.01],[0,0.3])
#plt.contour(f_sub,mass,np.transpose(pdf))
plt.contour(f_sub,mass,np.transpose(pdf_com))
plt.xlabel('log10(fsub)')
plt.ylabel('log10(half-mode mass)')
#plt.title(lens)
plt.title('5 lenses joint contour')
plt.ylim(6,9)
#plt.twinx()
#plt.legend(scatterpoints=1)
plt.colorbar()
#plt.show()
#plt.savefig(path+'lens_info/'+lens_list+'_wdm_01.png')
plt.savefig(path+'lens_info/5lens_wdm.png')
