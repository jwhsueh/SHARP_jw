import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

lens='mock3'
path = '/Volumes/sting_1/subs/'+lens+'/protect/'


### fsub=0.01
#t1 = np.loadtxt(path+lens+'_abc0100_chi2_com.txt')
t1 = np.loadtxt(path+lens+'_abc0100_0chi2.txt')
t0100 = t1

## fsub=0.005
#t1 = np.loadtxt(path+lens+'_abc00500_chi2_com.txt')
t1 = np.loadtxt(path+lens+'_abc00500_0chi2.txt')
#t2 = np.loadtxt(path+'B1422_abc00500_chi2_com.txt')
#t3 = np.loadtxt(path+'B1422_abc00500_chi2_com2.txt')
#t00500 = np.append(t1,t2)
#t00500 = np.append(t00500,t3)
t00500=t1
'''
## fsub=0.002
t1 = np.loadtxt(path+lens+'_abc00200_chi2_com.txt')
#t1 = np.loadtxt(path+lens+'_abc00200_0chi2.txt')
t00200 = t1+1

## fsub=0.001
#t1 = np.loadtxt(path+lens+'_abc00100_1chi2.txt')
t1 = np.loadtxt(path+lens+'_abc00100_chi2_com.txt')
t00100 = t1+1
'''

## fsub=0.003
#t1 = np.loadtxt(path+'mock2_abc00300_0chi2.txt')
#t00300 = t1



## fsub=0.0005
#t1 = np.loadtxt(path+'mock2_abc000500_0chi2.txt')
#t000500 = t1


## fsub=0.02
#t1 = np.loadtxt(path+lens+'_abc0200_chi2_com.txt')
t1 = np.loadtxt(path+lens+'_abc0200_0chi2.txt')
#t3 = np.loadtxt(path+'B1422_abc0200_chi2_com2.txt')
t0200 = t1

## fsub=0.03
#t1 = np.loadtxt(path+lens+'_abc0300_chi2_com.txt')
t1 = np.loadtxt(path+lens+'_abc0300_0chi2.txt')
#t1 = np.loadtxt(path+lens+'_abc0200_0chi2.txt')
#t3 = np.loadtxt(path+'B1422_abc0200_chi2_com2.txt')
t0300 = t1

## fsub=0.05
t1 = np.loadtxt(path+lens+'_abc0500_chi2_com.txt')
#t2 = np.loadtxt(path+'B1422_abc0500_1chi2.txt')
#t3 = np.loadtxt(path+'B1422_abc0500_2chi2.txt')
t0500 = t1
#t0500 = np.append(t0500,t3)
'''
## fsub=0.1
t1 = np.loadtxt(path+lens+'_abc1000_chi2_com.txt')
#t1 = np.loadtxt(path+lens+'_abc1000_0chi2.txt')
t1000 = t1
#t1000 = np.append(t1000,t3)
#t1000 = np.append(t1000,t4)
'''
## fsub=0.2
#t1 = np.loadtxt(path+'B1422_abc2000_0chi2.txt')
#t2 = np.loadtxt(path+'B1422_abc2000_1chi2.txt')
#t2000 = np.append(t1,t2)
###

'''
th = 0.01

n0100 =int(len(t0100)*th)
n00500 =int(len(t00500)*th)
n0200 =int(len(t0200)*th)
n0500 =int(len(t0500)*th)

th0100 = np.sort(t0100)[n0100]
th00500 = np.sort(t00500)[n00500]
th0200 = np.sort(t0200)[n0200]
th0500 = np.sort(t0500)[n0500]

t0100 = t0100[t0100<th0100]
t00500 = t00500[t00500<th00500]
t0200 = t0200[t0200<th0200]
t0500 = t0500[t0500<th0500]
'''
###
sn = np.array([1000,5000,1e4,1.5e4,3e4,5e4])
#sn = np.array([30,100,300,500,1000])
#sn = np.array([500,1000,5000,3.5e4])

#t = t0100
select = sn[3]
N=100
#p0100,p1000,p0200,p0500,p00500,p00200,p00100 = np.zeros(N),np.zeros(N),np.zeros(N),np.zeros(N),np.zeros(N),np.zeros(N),np.zeros(N)
#p0100,p0200,p00500,p00200,p00100,p0500,p1000 = np.zeros(N),np.zeros(N),np.zeros(N),np.zeros(N),np.zeros(N),np.zeros(N),np.zeros(N)
p0100,p0200,p0300,p00500,p0500 = np.zeros(N),np.zeros(N),np.zeros(N),np.zeros(N),np.zeros(N)
#p00500,p00300,p00200,p00100,p000500 = np.zeros(N),np.zeros(N),np.zeros(N),np.zeros(N),np.zeros(N)

for i in range(N):
	pick = np.random.choice(t0100,select,replace=False)
	p0100[i] = np.sum(np.exp(-pick/2.))
	pick = np.random.choice(t0200,select,replace=False)
	p0200[i] = np.sum(np.exp(-pick/2.))
	pick = np.random.choice(t0500,select,replace=False)
	p0300[i] = np.sum(np.exp(-pick/2.))
	pick = np.random.choice(t0300,select,replace=False)
	p0500[i] = np.sum(np.exp(-pick/2.))
	#pick = np.random.choice(t1000,select,replace=False)
	#p1000[i] = np.sum(np.exp(-pick/2.))
	pick = np.random.choice(t00500,select,replace=False)
	p00500[i] = np.sum(np.exp(-pick/2.))
	#pick = np.random.choice(t00200,select,replace=False)
	#p00200[i] = np.sum(np.exp(-pick/2.))
	#pick = np.random.choice(t00300,select,replace=False)
	#p00300[i] = np.sum(np.exp(-pick/2.))
	#pick = np.random.choice(t00100,select,replace=False)
	#p00100[i] = np.sum(np.exp(-pick/2.))
	#pick = np.random.choice(t000500,select,replace=False)
	#p000500[i] = np.sum(np.exp(-pick/2.))

'''
#f_series = np.array(np.log10([0.005,0.003,0.002,0.001]))
f_series = np.array(np.log10([0.005,0.01,0.02,0.03,0.05]))
#f_series = np.array(np.log10([0.005,0.01,0.02,0.05,0.1]))
#f_series = np.array([0.001,0.002,0.005,0.01,0.02,0.05,0.1])
#p_series = np.array([np.average(p00500),np.average(p0100),np.average(p0200),np.average(p0500),np.average(p1000)])
#p_series = np.array([np.average(p00500),np.average(p00300),np.average(p00200),np.average(p00100)])
p_series = np.array([np.average(p00500),np.average(p0100),np.average(p0200),np.average(p0300),np.average(p0500)])
#p_series = np.array([np.average(p00100),np.average(p00200),np.average(p00500),np.average(p0100),np.average(p0200),np.average(p0500),np.average(p1000)])
#y_err = [np.std(p00100),np.std(p00200),np.std(p00500),np.std(p0100),np.std(p0200),np.std(p0500),np.std(p1000)]
#y_err = [np.std(p00500),np.std(p0100),np.std(p0200),np.std(p0500),np.std(p1000)]
#y_err = [np.std(p00500),np.std(p00300),np.std(p00200),np.std(p00100)]
y_err = [np.std(p00500),np.std(p0100),np.std(p0200),np.std(p0300),np.std(p0500)]
plt.errorbar(f_series,p_series,yerr = y_err,marker='o')

#p_series = np.array([np.average(p0100),np.average(p00500),np.average(p00300),np.average(p00200),np.average(p00100)])
#y_err = [np.std(p0100),np.std(p00500),np.std(p00300),np.std(p00200),np.std(p00100)]
#plt.errorbar(f_series,p_series,yerr = y_err,marker='o')
#p_series = p_series/np.max(p_series)

#plt.semilogx(f_series,p_series,marker='.',linestyle = 'none')
#plt.plot([np.log10(0.00314),np.log10(0.00314)],[0,2],c='r',label='real f_sub = 0.314%')
plt.plot([np.log10(0.0229),np.log10(0.0229)],[2,3],c='r',label='real f_sub = 2.29%')
'''

def gauss(x,a,x0,sig):
	return a*np.exp(-(x-x0)**2/(2*sig**2))
'''
mean,sig = np.log10(0.00314),np.log10(0.0001)

popt,pcov = curve_fit(gauss,f_series,p_series,p0=[1,mean,sig],sigma = y_err)

x = np.linspace(min(f_series), max(f_series), 100)
y = gauss(x, *popt)
plt.plot(x, y)
'''
'''
#plt.xlim(0.0,0.25)
#plt.xlim(-3.1,-2.2)
plt.ylim(2,3)

plt.xlabel('log(f_sub)')
plt.ylabel('sum of likelihood')
plt.title('mock3')
plt.legend(loc=2)
#plt.show()

'''
plt.title('mock3, Sum of likelihood of 15000 realizations drawn from entire pool')

plt.semilogy(np.arange(N),p0500,c='r',label='fsub=5%')
plt.semilogy(np.arange(N),p0100,c='b',label='fsub=1%')
plt.semilogy(np.arange(N),p0200,c='k',label='fsub=2%')
plt.semilogy(np.arange(N),p0300,c='g',label='fsub=3%')
#plt.semilogy(np.arange(N),p1000,c='g',label='fsub=10%')
plt.semilogy(np.arange(N),p00500,c='r',linestyle='dotted',label='fsub=0.5%')
#plt.semilogy(np.arange(N),p00200,c='b',linestyle='dotted',label='fsub=0.2%')
#plt.semilogy(np.arange(N),p00300,c='b',label='fsub=0.3%')
#plt.semilogy(np.arange(N),p00100,c='k',linestyle='dotted',label='fsub=0.1%')
plt.ylabel('likelihood')
plt.legend()

#plt.show()

plt.savefig(path+'mock3_sol.png')

