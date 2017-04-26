import numpy as np
import matplotlib.pyplot as plt

path='/Volumes/sting_1/data/'

list_file = path+'snap99_elp2_Rcusp.txt'
#list_file_sub = path+'snap99_elp2_sub_Rcusp.txt'
list_file_sub = path+'snap99_tri2_Rcusp.txt'  # disc
list_file_in = path+'snap99_edg_size.txt'
size_tab = np.loadtxt(list_file_in)

file_list = np.genfromtxt(list_file,dtype='str')
file_list_sub = np.genfromtxt(list_file_sub,dtype='str')
theta = np.loadtxt(list_file_in)[:,2]

Rfold,Rcusp,phi0,phi1 = np.empty(1),np.empty(1),np.empty(1),np.empty(1)
for file_one in file_list:
	#print file_one
	table = np.loadtxt(path+'Rcusp_f/'+file_one+'_64_Rcusp_ga.txt')
	#mask = np.loadtxt(path+'cm_test/snap99_'+file_one+'_hfmask.txt').astype(bool)
	#mask = mask[:(table[:,0]).size]
	#Rfold,Rcusp,phi0,phi1 = np.append(Rfold,table[:,0][mask]),np.append(Rcusp,table[:,1][mask]),np.append(phi0,table[:,2][mask]),np.append(phi1,table[:,3][mask])
	Rfold,Rcusp,phi0,phi1 = np.append(Rfold,table[:,0]),np.append(Rcusp,table[:,1]),np.append(phi0,table[:,2]),np.append(phi1,table[:,3])


Rfold_s,Rcusp_s,phi0_s,phi1_s  = np.empty(1),np.empty(1),np.empty(1),np.empty(1)


for file_one in file_list_sub:
	#print file_one
	table = np.loadtxt(path+'Rcusp_f/'+file_one+'_64_Rcusp_ga.txt')
	#mask = size_tab[:,3]>0.25
	#mask = mask[:(table[:,0]).size]
	#Rfold_s,Rcusp_s,phi0_s,phi1_s = np.append(Rfold_s,table[:,0][mask]),np.append(Rcusp_s,table[:,1][mask]),np.append(phi0_s,table[:,2][mask]),np.append(phi1_s,table[:,3][mask])
	Rfold_s,Rcusp_s,phi0_s,phi1_s = np.append(Rfold_s,table[:,0]),np.append(Rcusp_s,table[:,1]),np.append(phi0_s,table[:,2]),np.append(phi1_s,table[:,3])


## edge-on selection
'''
for i in range(file_list_sub.size):
	file_one = file_list_sub[i]
	ina = theta[i]

	if np.logical_and(ina>=80,ina<=100):
		table = np.loadtxt(path+'Rcusp_c/'+file_one+'_64_Rcusp_ga.txt')
		Rfold_s,Rcusp_s,phi0_s,phi1_s = np.append(Rfold_s,table[:,0]),np.append(Rcusp_s,table[:,1]),np.append(phi0_s,table[:,2]),np.append(phi1_s,table[:,3])
'''

mask = np.abs(Rfold)<0.9
Rfold,Rcusp,phi0,phi1=Rfold[mask],Rcusp[mask],phi0[mask],phi1[mask]
mask = np.abs(Rcusp)<0.9
Rfold,Rcusp,phi0,phi1=Rfold[mask],Rcusp[mask],phi0[mask],phi1[mask]

mask = np.abs(Rfold_s)<0.7
Rfold_s,Rcusp_s,phi0_s,phi1_s=Rfold_s[mask],Rcusp_s[mask],phi0_s[mask],phi1_s[mask]
mask = np.abs(Rcusp_s)<0.7
Rfold_s,Rcusp_s,phi0_s,phi1_s=Rfold_s[mask],Rcusp_s[mask],phi0_s[mask],phi1_s[mask]


#Rcusp_s=Rcusp[np.abs(Rcusp)<1.0]

#Rfold=Rfold[Rfold<0.5]
#Rfold_s=Rfold_s[Rfold_s<0.5]

print Rfold.size, Rfold_s.size,phi0.size,phi1_s.size
#### 

## ---- probability contour
'''
## cusp
edge = np.linspace(30,120,10)
print edge

prob_table = np.zeros((5,edge.size-1))
prob_list = np.array([0.5,0.2,0.1,0.05,0.01])
prob_sig = np.zeros((5,edge.size-1))

for j in range(edge.size-1):
	mask = np.logical_and(phi0_s>=edge[j],phi0_s<edge[j+1])
	#mask = np.logical_and(phi0>=edge[j],phi0<edge[j+1])

	#phi0_b = phi0[mask]
	Rcusp_b = np.abs(Rcusp_s[mask])
	#Rcusp_b = np.abs(Rcusp[mask])
	print Rcusp_b.size
	sig = np.sqrt(Rcusp_b.size)

	Rcusp_b = np.sort(Rcusp_b)
	idx = np.round(Rcusp_b.size*prob_list).astype(int)

	prob_table[:,j] = Rcusp_b[-idx]
	prob_sig[:,j] = Rcusp_b[-idx]/sig
	print Rcusp_b[-idx]

## fold
'''
edge2 = np.linspace(0,50,11)
print edge2

prob_table = np.zeros((5,edge2.size-1))
prob_sig = np.zeros((5,edge2.size-1))
#prob_sig = np.zeros(edge2.size-1)
prob_list = np.array([0.5,0.2,0.1,0.05,0.01])

for j in range(edge2.size-1):
	mask = np.logical_and(phi1_s>=edge2[j],phi1_s<edge2[j+1])
	#mask = np.logical_and(phi1/57.3>=edge2[j],phi1/57.3<edge2[j+1])

	Rfold_b = np.abs(Rfold_s[mask])
	#Rfold_b = np.abs(Rfold[mask])
	print Rfold_b.size, np.sqrt(Rfold_b.size)
	sig = np.sqrt(Rfold_b.size)
	#sig_list = prob_list+sig
	print np.average(Rfold_b),sig*np.average(Rfold_b)
	#print np.sqrt(Rfold_b.size)/Rfold_b.size

	Rfold_b = np.sort(Rfold_b)
	idx = np.round(Rfold_b.size*prob_list).astype(int)
	#idx_p = np.round(Rfold_b.size*sig_list).astype(int)
	#print idx_p

	prob_table[:,j] = Rfold_b[-idx]
	prob_sig[:,j] = Rfold_b[-idx]/sig
	#prob_sig[:,j] = Rfold_b[-idx]-Rfold_b[-idx_p]
	print Rfold_b[-idx]
	#print Rfold_b[-idx]-Rfold_b[-idx_p]


####
## plot data points
plt.scatter(phi1_s,np.abs(Rfold_s),marker='x',color='r',alpha=0.3)
#plt.scatter(phi0_s,np.abs(Rcusp_s),marker='x',color='r',alpha=0.3)

####
'''
# cusp
plt.plot(edge[0:-1]+5,prob_table[4,0:],color='k',linestyle='-.',label='1%')
plt.plot(edge[0:-1]+5,prob_table[3,:],color='g',label='5%')
plt.plot(edge[0:-1]+5,prob_table[2,:],color='r',label='10%')
plt.plot(edge[0:-1]+5,prob_table[1,:],color='b',label='20%')
plt.plot(edge[0:-1]+5,prob_table[0,:],color='k',label='50%')
plt.fill_between(edge[0:-1]+5,prob_table[4,0:]-prob_sig[4,0:],prob_table[4,0:]+prob_sig[4,0:],color='k',alpha=0.1)
plt.fill_between(edge[0:-1]+5,prob_table[2,:]-prob_sig[2,:],prob_table[2,:]+prob_sig[2,:],color='r',alpha=0.1)
plt.fill_between(edge[0:-1]+5,prob_table[3,:]-prob_sig[3,:],prob_table[3,:]+prob_sig[3,:],color='g',alpha=0.3)
plt.fill_between(edge[0:-1]+5,prob_table[0,:]-prob_sig[0,:],prob_table[0,:]+prob_sig[0,:],color='k',alpha=0.3)
plt.fill_between(edge[0:-1]+5,prob_table[1,:]-prob_sig[1,:],prob_table[1,:]+prob_sig[1,:],color='b',alpha=0.3)

'''
# fold
plt.plot(edge2[0:-1]+2.5,prob_table[4,:],color='k',linestyle='-.',label='1%')
plt.plot(edge2[0:-1]+2.5,prob_table[3,:],color='g',label='5%')
plt.plot(edge2[0:-1]+2.5,prob_table[2,:],color='r',label='10%')
plt.plot(edge2[0:-1]+2.5,prob_table[1,:],color='b',label='20%')
plt.plot(edge2[0:-1]+2.5,prob_table[0,:],color='k',label='50%')
plt.fill_between(edge2[0:-1]+2.5,prob_table[4,:]-prob_sig[4,:],prob_table[4,:]+prob_sig[4,:],color='k',alpha=0.1)
plt.fill_between(edge2[0:-1]+2.5,prob_table[2,:]-prob_sig[2,:],prob_table[2,:]+prob_sig[2,:],color='r',alpha=0.1)
plt.fill_between(edge2[0:-1]+2.5,prob_table[3,:]-prob_sig[3,:],prob_table[3,:]+prob_sig[3,:],color='g',alpha=0.3)
plt.fill_between(edge2[0:-1]+2.5,prob_table[0,:]-prob_sig[0,:],prob_table[0,:]+prob_sig[0,:],color='k',alpha=0.3)
plt.fill_between(edge2[0:-1]+2.5,prob_table[1,:]-prob_sig[1,:],prob_table[1,:]+prob_sig[1,:],color='b',alpha=0.3)


## flux ratio data point
tab = np.loadtxt('../../data/flux_ratio.txt')
tab[:,1],tab[:,3] = np.abs(tab[:,1]),np.abs(tab[:,3])
'''
#cusp
plt.scatter(tab[:,0],tab[:,1],color='k',marker='^',s=100)
plt.text(tab[1,0]-4,tab[1,1]-0.05,'MG0414',color='k',fontsize=14)
plt.text(tab[2,0]+1,tab[2,1]-0.01,'B1422',color='k',fontsize=14)
plt.text(tab[3,0]-4,tab[3,1]+0.02,'B1608',color='k',fontsize=14)
plt.plot(43,tab[4,1],marker='o',mec='k',mfc='none',ms=6,mew=1.3)
plt.text(42,tab[4,1]-0.04,'B2045',color='k',fontsize=14)
plt.errorbar(42.5,tab[4,1]+0.002,xerr=2.5,xuplims=True,color='k',elinewidth=2,capsize=4)
'''
#fold
plt.scatter(tab[:,2],tab[:,3],color='k',marker='^',s=100)
plt.text(tab[0,2]-2,tab[0,3]+0.02,'B0128',color='k',fontsize=14)
plt.text(tab[1,2]+0.5,tab[1,3],'MG0414',color='k',fontsize=14)
plt.text(tab[4,2]-2,tab[4,3]+0.03,'B2045',color='k',fontsize=14)

tab = np.loadtxt('../../data/flux_ratio_disk.txt')
tab[:,1],tab[:,3] = np.abs(tab[:,1]),np.abs(tab[:,3])
'''
#cusp
plt.scatter(tab[:,0],tab[:,1],color='b',marker='d',s=100)
plt.text(tab[0,0]+1,tab[0,1]-0.01,'B0712',color='b',fontsize=14)
plt.text(tab[1,0]-9.5,tab[1,1]-0.01,'B1555',color='b',fontsize=14)
'''
#fold
plt.scatter(tab[:,2],np.abs(tab[:,3]),color='b',marker='d',s=100)
plt.text(tab[0,2]-5,tab[0,3]-0.01,'B0712',color='b',fontsize=14)
plt.text(tab[1,2]+0.5,tab[1,3]+0.01,'B1555',color='b',fontsize=14)

##

#plt.text(42,0.85,'Elliptical lens',fontsize=14)
plt.text(6,0.85,'Disc lens',fontsize=14)
#plt.text(42,0.8,'$P$ $(>|R_{cusp}|)$',fontsize=14)
plt.text(6,0.8,'$P$ $(>|R_{fold}|)$',fontsize=14)

#plt.title('edge-on disks')
#plt.title('elliptical gauss fit')
#plt.ylabel('|$R_{cusp}$|')
plt.ylabel('|$R_{fold}$|')
#plt.ylabel('number counts')
#plt.xlabel('$\Delta \phi$')
plt.xlabel('$\phi_1$')
#plt.xlabel('theta1/theta_E')
#plt.legend(loc=1)
#plt.xlim(40,110)
plt.xlim(5,40)
plt.ylim(0,0.9)
#plt.gca().set_aspect(70/0.9)
plt.gca().set_aspect(35/0.9)
#plt.show()
#plt.savefig('../../data/glamer/glamer_tri_fold_pd.png',bbox_inches='tight')

