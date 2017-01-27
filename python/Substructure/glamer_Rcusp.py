import numpy as np
import matplotlib.pyplot as plt

path='/Volumes/sting_1/data/'

list_file = path+'snap99_elp2_sub_Rcusp.txt'
#list_file_sub = path+'snap99_elp2_sub_Rcusp.txt'
list_file_sub = path+'snap99_tri2_sub_Rcusp.txt'  # disc

file_list = np.genfromtxt(list_file,dtype='str')
file_list_sub = np.genfromtxt(list_file_sub,dtype='str')

Rfold,Rcusp = np.empty(1),np.empty(1)
for file_one in file_list:
	table = np.loadtxt(path+'Rcusp/'+file_one+'.txt')
	Rfold,Rcusp = np.append(Rfold,table[:,0]),np.append(Rcusp,table[:,1])

Rfold_s,Rcusp_s = np.empty(1),np.empty(1)
for file_one in file_list_sub:
	table = np.loadtxt(path+'Rcusp/'+file_one+'.txt')
	Rfold_s,Rcusp_s = np.append(Rfold_s,table[:,0]),np.append(Rcusp_s,table[:,1])

Rfold=Rfold[np.abs(Rfold)<1.0]
Rfold_s=Rfold_s[np.abs(Rfold_s)<1.0]
#Rcusp_s=Rcusp[np.abs(Rcusp)<1.0]

#Rfold=Rfold[Rfold<0.5]
#Rfold_s=Rfold_s[Rfold_s<0.5]

print Rfold.size, Rfold_s.size
#### 
'''
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
'''
w_Rfold = np.ones_like(Rfold)/len(Rfold)
w_Rfold_s = np.ones_like(Rfold_s)/len(Rfold_s)
w_Rcusp = np.ones_like(Rcusp)/len(Rcusp)
w_Rcusp_s = np.ones_like(Rcusp_s)/len(Rcusp_s)

#plt.hist(np.abs(Rcusp),alpha=0.5,color='b',label='non disc',weights=w_Rcusp,bins=np.linspace(0,0.6,30))
#plt.hist(np.abs(Rcusp_s),alpha=0.5,color='r',label='disc',weights=w_Rcusp_s,bins=np.linspace(0,0.6,30))
#plt.xlabel('|Rcusp|')
plt.title('w/ sub')

plt.hist(np.abs(Rfold),alpha=0.5,color='b',label='non disk',weights=w_Rfold,bins=np.linspace(0,0.4,40))
plt.hist(np.abs(Rfold_s),alpha=0.5,color='r',label='disk',weights=w_Rfold_s,bins=np.linspace(0,0.4,40))
plt.xlabel('|Rfold|')
plt.legend()
#plt.show()
plt.savefig('../../data/glamer/glamer_Rfold_sub_half.png')

