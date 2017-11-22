import numpy as np
import matplotlib.pyplot as plt

path = '/Volumes/sting_1/subs/B1422/protect/'


### fsub=0.01
t1 = np.loadtxt(path+'B1422_abc0100_0chi2.txt')
t2 = np.loadtxt(path+'B1422_abc0100_1chi2.txt')
t3 = np.loadtxt(path+'B1422_abc0100_chi2_com.txt')
t0100 = np.append(t1,t2)
t0100 = np.append(t0100,t3)

## fsub=0.005
t1 = np.loadtxt(path+'B1422_abc00500_0chi2.txt')
t2 = np.loadtxt(path+'B1422_abc00500_chi2_com.txt')
t3 = np.loadtxt(path+'B1422_abc00500_chi2_com2.txt')
t00500 = np.append(t1,t2)
t00500 = np.append(t00500,t3)

## fsub=0.02
t2 = np.loadtxt(path+'B1422_abc0200_chi2_com.txt')
t3 = np.loadtxt(path+'B1422_abc0200_chi2_com2.txt')
t0200 = np.append(t2,t3)

## fsub=0.05
t1 = np.loadtxt(path+'B1422_abc0500_0chi2.txt')
t2 = np.loadtxt(path+'B1422_abc0500_1chi2.txt')
t3 = np.loadtxt(path+'B1422_abc0500_2chi2.txt')
t0500 = np.append(t1,t2)
t0500 = np.append(t0500,t3)

###
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

###
#sn = np.array([1000,5000,1e4,5e4])
sn = np.array([30,100,300,500,1000])

#t = t0100
select = sn[3]
N=100
p0100,p00500,p0200,p0500 = np.zeros(N),np.zeros(N),np.zeros(N),np.zeros(N)

for i in range(N):
	pick = np.random.choice(t0100,select,replace=False)
	p0100[i] = np.sum(np.exp(-pick/2.))
	pick = np.random.choice(t00500,select,replace=False)
	p00500[i] = np.sum(np.exp(-pick/2.))
	pick = np.random.choice(t0200,select,replace=False)
	p0200[i] = np.sum(np.exp(-pick/2.))
	pick = np.random.choice(t0500,select,replace=False)
	p0500[i] = np.sum(np.exp(-pick/2.))

plt.title('Sum of likelihood of 500 realizations drawn from top 1% pool')

plt.semilogy(np.arange(N),p00500,c='r',label='fsub=0.5%')
plt.semilogy(np.arange(N),p0100,c='b',label='fsub=1%')
plt.semilogy(np.arange(N),p0200,c='k',label='fsub=2%')
plt.semilogy(np.arange(N),p0500,c='g',label='fsub=5%')
plt.ylabel('likelihood')
plt.legend()
plt.gca().set_axis(0.5/0.7)
plt.show()
#plt.savefig(path+'B1422_likelihood_500p.png')

