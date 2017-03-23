import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.ndimage.filters import gaussian_filter
from scipy import stats

path='/Volumes/sting_1/data/'

list_file = path+'snap99_elp2_Rcusp.txt'
#list_file_sub = path+'snap99_elp2_sub_Rcusp.txt'
list_file_sub = path+'snap99_tri2_Rcusp.txt'  # disc
list_file_in = path+'snap99_tri2_size.txt'

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

'''
for file_one in file_list_sub:
	#print file_one
	table = np.loadtxt(path+'Rcusp_c/'+file_one+'_64_Rcusp_ga.txt')
	#mask = np.loadtxt(path+'cm_test/snap99_'+file_one+'_hfmask.txt').astype(bool)
	#mask = mask[:(table[:,0]).size]
	#Rfold_s,Rcusp_s,phi0_s,phi1_s = np.append(Rfold_s,table[:,0][mask]),np.append(Rcusp_s,table[:,1][mask]),np.append(phi0_s,table[:,2][mask]),np.append(phi1_s,table[:,3][mask])
	Rfold_s,Rcusp_s,phi0_s,phi1_s = np.append(Rfold_s,table[:,0]),np.append(Rcusp_s,table[:,1]),np.append(phi0_s,table[:,2]),np.append(phi1_s,table[:,3])
'''
for i in range(file_list_sub.size):
	file_one = file_list_sub[i]
	ina = theta[i]

	if np.logical_and(ina>=80,ina<=100):
		table = np.loadtxt(path+'Rcusp_c/'+file_one+'_64_Rcusp_ga.txt')
		Rfold_s,Rcusp_s,phi0_s,phi1_s = np.append(Rfold_s,table[:,0]),np.append(Rcusp_s,table[:,1]),np.append(phi0_s,table[:,2]),np.append(phi1_s,table[:,3])


mask = np.abs(Rfold)<0.9
Rfold,Rcusp,phi0,phi1=Rfold[mask],Rcusp[mask],phi0[mask],phi1[mask]
mask = np.abs(Rcusp)<0.9
Rfold,Rcusp,phi0,phi1=Rfold[mask],Rcusp[mask],phi0[mask],phi1[mask]

mask = np.abs(Rfold_s)<0.9
Rfold_s,Rcusp_s,phi0_s,phi1_s=Rfold_s[mask],Rcusp_s[mask],phi0_s[mask],phi1_s[mask]
mask = np.abs(Rcusp_s)<0.9
Rfold_s,Rcusp_s,phi0_s,phi1_s=Rfold_s[mask],Rcusp_s[mask],phi0_s[mask],phi1_s[mask]

## inclination mask
mask = np.abs(Rfold_s)<0.9
Rfold_s,Rcusp_s,phi0_s,phi1_s=Rfold_s[mask],Rcusp_s[mask],phi0_s[mask],phi1_s[mask]
mask = np.abs(Rcusp_s)<0.9
Rfold_s,Rcusp_s,phi0_s,phi1_s=Rfold_s[mask],Rcusp_s[mask],phi0_s[mask],phi1_s[mask]

#Rcusp_s=Rcusp[np.abs(Rcusp)<1.0]

#Rfold=Rfold[Rfold<0.5]
#Rfold_s=Rfold_s[Rfold_s<0.5]

print Rfold.size, Rfold_s.size,phi0.size,phi1_s.size
#### 

## ---- probability contour
## cusp
edge = np.linspace(50,150,11)

prob_table = np.zeros((5,edge.size-1))
prob_list = np.array([0.5,0.2,0.1,0.05,0.01])

for j in range(edge.size-1):
	#mask = np.logical_and(phi0_s>=edge[j],phi0_s<edge[j+1])
	mask = np.logical_and(phi0>=edge[j],phi0<edge[j+1])

	#phi0_b = phi0[mask]
	#Rcusp_b = np.abs(Rcusp_s[mask])
	Rcusp_b = np.abs(Rcusp[mask])
	print Rcusp_b.size

	Rcusp_b = np.sort(Rcusp_b)
	idx = np.round(Rcusp_b.size*prob_list).astype(int)

	prob_table[:,j] = Rcusp_b[-idx]
	print Rcusp_b[-idx]

## fold

edge2 = np.linspace(0.2,0.8,7)
print edge2

prob_table = np.zeros((5,edge2.size-1))
prob_list = np.array([0.5,0.2,0.1,0.05,0.01])

for j in range(edge2.size-1):
	mask = np.logical_and(phi1_s/57.3>=edge2[j],phi1_s/57.3<edge2[j+1])
	#mask = np.logical_and(phi1/57.3>=edge2[j],phi1/57.3<edge2[j+1])

	Rfold_b = np.abs(Rfold_s[mask])
	#Rfold_b = np.abs(Rfold[mask])
	print Rfold_b.size

	Rfold_b = np.sort(Rfold_b)
	idx = np.round(Rfold_b.size*prob_list).astype(int)

	prob_table[:,j] = Rfold_b[-idx]
	print Rfold_b[-idx]


####
w_Rfold = np.ones_like(Rfold)/len(Rfold)
w_Rfold_s = np.ones_like(Rfold_s)/len(Rfold_s)
w_Rcusp = np.ones_like(Rcusp)/len(Rcusp)
w_Rcusp_s = np.ones_like(Rcusp_s)/len(Rcusp_s)


#plt.hist(np.abs(Rcusp),alpha=0.5,color='b',label='non disc',weights=w_Rcusp,bins=np.linspace(0,0.6,30))
#plt.hist(np.abs(Rcusp_s),alpha=0.5,color='r',label='disc',weights=w_Rcusp_s,bins=np.linspace(0,0.6,30))
#plt.xlabel('|Rcusp|')
#plt.title('w/ sub')

#plt.hist(np.abs(Rfold),alpha=0.5,color='b',label='non disk',weights=w_Rfold,bins=np.linspace(0,0.4,40))
#plt.hist(np.abs(Rfold_s),alpha=0.5,color='r',label='disk',weights=w_Rfold_s,bins=np.linspace(0,0.4,40))
#plt.xlabel('|Rfold|')

#plt.hist(phi0_s,bins=np.linspace(0,180,40))



#H,xbin,ybin = np.histogram2d(phi0,np.abs(Rcusp),bins=(np.linspace(0,180,40),np.linspace(0,0.5,50)))
H,xbin,ybin = np.histogram2d(phi1,(Rfold),bins=(np.linspace(0,60,40),np.linspace(0,0.5,50)))
H=H.T
#H = gaussian_filter(H,2.0)
#H = np.abs(H -1)
plt.imshow(H,interpolation='nearest',origin='low',extent=[xbin[0],xbin[-1],ybin[0],ybin[-1]])
#fig = plt.contour(xbin[:-1],ybin[:-1],H,[1,4,7,10,13])
#fig = plt.contour(xbin[:-1],ybin[:-1],H,[1,2,3,4,5])
#fig = plt.contour(xbin[:-1],ybin[:-1],H)
#plt.clabel(fig)
plt.gca().set_aspect(60/0.4)
#plt.gca().set_aspect(180/0.5)

'''
## try w/ KDE
values = np.vstack([phi0,Rcusp])
X, Y = np.mgrid[xbin[0]:xbin[-1]:40j,ybin[0]:ybin[-1]:50j]
positions = np.vstack([X.ravel(), Y.ravel()])
kernel = stats.gaussian_kde(values)
Z = np.reshape(kernel(positions).T, X.shape)
print np.max(Z), Z.shape
#plt.imshow(np.rot90(Z)/np.max(Z),interpolation='nearest',extent=[xbin[0],xbin[-1],ybin[0],ybin[-1]])
fig = plt.contour(X,Y,Z/np.max(Z))
plt.clabel(fig)
'''


'''
plt.plot(edge[0:-1]+5,prob_table[4,:],color='k',linestyle='-.',label='1%')
plt.plot(edge[0:-1]+5,prob_table[3,:],color='g',label='5%')
plt.plot(edge[0:-1]+5,prob_table[2,:],color='r',label='10%')
plt.plot(edge[0:-1]+5,prob_table[1,:],color='b',label='20%')
plt.plot(edge[0:-1]+5,prob_table[0,:],color='k',label='50%')


plt.plot(edge2[0:-1]+5/57.3,prob_table[4,:],color='k',linestyle='-.',label='1%')
plt.plot(edge2[0:-1]+5/57.3,prob_table[3,:],color='g',label='5%')
plt.plot(edge2[0:-1]+5/57.3,prob_table[2,:],color='r',label='10%')
plt.plot(edge2[0:-1]+5/57.3,prob_table[1,:],color='b',label='20%')
plt.plot(edge2[0:-1]+5/57.3,prob_table[0,:],color='k',label='50%')
'''

#B1555 & B0712
#plt.scatter([0.365,0.243],[0.235,0.085],marker='*',color='k',s=100)

#plt.title('disks half area (pt~750)')
plt.title('disks gauss fit ')
#plt.title('elliptical gauss fit')
#plt.ylabel('|Rcusp|')
plt.ylabel('|Rfold|')
#plt.ylabel('number counts')
#plt.xlabel('delta phi')
#plt.xlabel('phi 1')
plt.xlabel('theta1/theta_E')
plt.legend(loc=1)
#plt.xlim(50,150)
#plt.xlim(10,50)
#plt.ylim(0,0.9)
#plt.colorbar()
plt.show()
#plt.savefig('../../data/glamer/glamer_tri_ga_fold_pd.png')

