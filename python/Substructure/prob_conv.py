import numpy as np
import matplotlib.pyplot as plt

chi2 = np.loadtxt('/Volumes/sting_1/subs/MG0414/result/MG0414_sub01_0chi2_copy2.txt')

chi2 = np.sort(chi2)
prob = np.exp(-1.0*chi2/2.0)

acc = np.zeros(len(prob))
idx = np.arange(len(prob))+1

a = 0.0
for i in range(len(prob)):
	a = a+prob[i]
	acc[i] = a

print chi2[30000],chi2[40000]
print len(prob)

'''
frac = np.linspace(0.01,0.5,100)
N = len(prob)*frac
sum_of = np.zeros(len(N))

for i in range(len(N)):
	p_list = np.random.choice(prob,N[i],replace='false')
	sum_of[i] = np.sum(p_list)/N[i]


#print len(sum_of), len(frac), len(N)
plt.plot(N,sum_of)
plt.xlabel('number of ray-tracing draw')
plt.ylabel('sum of likelihood/# ray-tracing')
plt.title('MG0414, N = 50,000')
plt.show()
#plt.savefig('MG0414_draw3.png')
'''

plt.plot(idx,np.log(acc))
plt.xlabel('log(number of ray-tracing)')
plt.ylabel('sum of likelihood')
plt.title('N = 10,000')
plt.show()
#plt.savefig('B1422_10k_acc.png')