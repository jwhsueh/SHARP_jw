import numpy as np
import matplotlib.pyplot as plt

path='/Volumes/sting_1/'

table_path=path+'snap99_263313_test/Rcusp_263313_p1sub_64.txt'
table=np.loadtxt(table_path)

Rfold_s,Rcusp_s = table[:,0],table[:,1]

table_path=path+'snap99_263313_test/Rcusp_263313_p2sub_64.txt'
table=np.loadtxt(table_path)

Rfold_s,Rcusp_s = np.append(Rfold_s,table[:,0]),np.append(Rcusp_s,table[:,1])

table_path=path+'snap99_269741_test/Rcusp_269741_p2sub_64.txt'
table=np.loadtxt(table_path)

Rfold_s,Rcusp_s = np.append(Rfold_s,table[:,0]),np.append(Rcusp_s,table[:,1])

table_path=path+'snap99_269779_test/Rcusp_269779_p3sub_64.txt'
table=np.loadtxt(table_path)

Rfold_s,Rcusp_s = np.append(Rfold_s,table[:,0]),np.append(Rcusp_s,table[:,1])

Rfold_s=Rfold_s[np.abs(Rfold_s)<0.1]
Rcusp_s=Rcusp_s[np.abs(Rcusp_s)<0.5]

#### 

table_path=path+'snap99_263313_test/Rcusp_263313_p1_64.txt'
table=np.loadtxt(table_path)

Rfold,Rcusp = table[:,0],table[:,1]

table_path=path+'snap99_263313_test/Rcusp_263313_p2_64.txt'
table=np.loadtxt(table_path)

Rfold,Rcusp = np.append(Rfold,table[:,0]),np.append(Rcusp,table[:,1])

table_path=path+'snap99_269741_test/269741_p2_64_Rcusp.txt'
table=np.loadtxt(table_path)

Rfold,Rcusp = np.append(Rfold,table[:,0]),np.append(Rcusp,table[:,1])

table_path=path+'snap99_269779_test/269779_p3_64_Rcusp.txt'
table=np.loadtxt(table_path)

Rfold,Rcusp = np.append(Rfold,table[:,0]),np.append(Rcusp,table[:,1])

Rfold=Rfold[np.abs(Rfold)<0.1]
Rcusp=Rcusp[np.abs(Rcusp)<0.5]

w_Rfold = np.ones_like(Rfold)/len(Rfold)
w_Rfold_s = np.ones_like(Rfold_s)/len(Rfold_s)
w_Rcusp = np.ones_like(Rcusp)/len(Rcusp)
w_Rcusp_s = np.ones_like(Rcusp_s)/len(Rcusp_s)

plt.hist(np.abs(Rcusp),alpha=0.5,color='b',label='no sub',weights=w_Rcusp)
plt.hist(np.abs(Rcusp_s),alpha=0.5,color='r',label='w/ sub',weights=w_Rcusp_s)
plt.xlabel('|Rcusp|')
plt.legend()
#plt.show()
plt.savefig('../../data/glamer/glamer_Rcusp_test_n.png')

