import numpy as np
import matplotlib.pyplot as plt

t1 = np.loadtxt('/Volumes/sting_1/subs/result/B1422_los00_0chi2.txt')
t2 = np.loadtxt('/Volumes/sting_1/subs/result/B1422_los00_2000chi2.txt')
t4 = np.loadtxt('/Volumes/sting_1/subs/result/B1422_los00_4000chi2.txt')
t6 = np.loadtxt('/Volumes/sting_1/subs/result/B1422_los00_6000chi2.txt')
t8 = np.loadtxt('/Volumes/sting_1/subs/result/B1422_los00_8000chi2.txt')

chi2 = np.append(t1,t2)
chi2 = np.append(chi2,t4)
chi2 = np.append(chi2,t6)
chi2 = np.append(chi2,t8)


s1 = np.loadtxt('/Volumes/sting_1/subs/result/B1422_sub01_0chi2.txt')
s2 = np.loadtxt('/Volumes/sting_1/subs/result/B1422_sub01_5000chi2.txt')
chi2_s = np.append(s1,s2)
'''
## randomly select chi2
n_line = 100
n_select = 500

for i in range(n_line):
	
	chi_sub = np.random.choice(chi2,n_select)
	chi_sub = chi_sub[chi_sub<50]

	hist_sub = np.histogram(chi_sub)
	node = hist_sub[1]

	plt.plot(node[:-1],hist_sub[0]/float(n_select),color='k',alpha=0.1)


plt.plot(node[:-1],hist_sub[0]/float(n_select),color='k',alpha=0.1,label='subset (500)')
'''

chi2_all = chi2[chi2<1e3]
'''
hist_all = np.histogram(chi2_all)
node = hist_all[1]
plt.plot(node[:-1],np.log10(hist_all[0]/float(10000)),color='r',label='only LOS')

chi2_sall = chi2_s[chi2_s<50]
hist_sall = np.histogram(chi2_sall)
node = hist_sall[1]
plt.plot(node[:-1],np.log10(hist_sall[0]/float(10000)),color='b',label='only subs, f_sub=1%')
'''

chi2_all.sort()
indx = np.arange(chi2_all.size).cumsum()+1.
indx /= indx[-1]
plt.plot(chi2_all,np.log10(indx),color='r',label = 'only LOS')

chi2_sall = chi2_s[chi2_s<1e3]
chi2_sall.sort()
indx = np.arange(chi2_sall.size).cumsum()+1.
indx /= indx[-1]
plt.plot(chi2_sall,np.log10(indx),label = 'only subs')

#plt.title('10000 realizations')
plt.xlabel('chi^2')
plt.ylabel('fraction')
plt.legend(loc=2)
plt.xlim(0,100)
#plt.show()
plt.savefig('sub_01_zoom.png')