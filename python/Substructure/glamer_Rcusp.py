import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

path='/Volumes/sting_1/data/'

list_file = path+'snap99_elp2_Rcusp.txt'
#list_file_sub = path+'snap99_elp2_sub_Rcusp.txt'
list_file_sub = path+'snap99_tri2_Rcusp.txt'  # disc

file_list = np.genfromtxt(list_file,dtype='str')
file_list_sub = np.genfromtxt(list_file_sub,dtype='str')

Rfold,Rcusp,phi0,phi1 = np.empty(1),np.empty(1),np.empty(1),np.empty(1)
for file_one in file_list:
	print file_one
	table = np.loadtxt(path+'Rcusp_c/'+file_one+'_64_Rcusp_s.txt')
	#mask = np.loadtxt(path+'cm_test/snap99_'+file_one+'_hfmask.txt').astype(bool)
	#mask = mask[:(table[:,0]).size]
	#Rfold,Rcusp,phi0,phi1 = np.append(Rfold,table[:,0][mask]),np.append(Rcusp,table[:,1][mask]),np.append(phi0,table[:,2][mask]),np.append(phi1,table[:,3][mask])
	Rfold,Rcusp,phi0,phi1 = np.append(Rfold,table[:,0]),np.append(Rcusp,table[:,1]),np.append(phi0,table[:,2]),np.append(phi1,table[:,3])

Rfold_s,Rcusp_s,phi0_s,phi1_s  = np.empty(1),np.empty(1),np.empty(1),np.empty(1)
for file_one in file_list_sub:
	print file_one
	table = np.loadtxt(path+'Rcusp_c/'+file_one+'_64_Rcusp_s.txt')
	#mask = np.loadtxt(path+'cm_test/snap99_'+file_one+'_hfmask.txt').astype(bool)
	#mask = mask[:(table[:,0]).size]
	#Rfold_s,Rcusp_s,phi0_s,phi1_s = np.append(Rfold_s,table[:,0][mask]),np.append(Rcusp_s,table[:,1][mask]),np.append(phi0_s,table[:,2][mask]),np.append(phi1_s,table[:,3][mask])
	Rfold_s,Rcusp_s,phi0_s,phi1_s = np.append(Rfold_s,table[:,0]),np.append(Rcusp_s,table[:,1]),np.append(phi0_s,table[:,2]),np.append(phi1_s,table[:,3])

mask = np.abs(Rfold)<0.4
Rfold,Rcusp,phi0,phi1=Rfold[mask],Rcusp[mask],phi0[mask],phi1[mask]
mask = np.abs(Rcusp)<0.4
Rfold,Rcusp,phi0,phi1=Rfold[mask],Rcusp[mask],phi0[mask],phi1[mask]

mask = np.abs(Rfold_s)<0.4
Rfold_s,Rcusp_s,phi0_s,phi1_s=Rfold_s[mask],Rcusp_s[mask],phi0_s[mask],phi1_s[mask]
mask = np.abs(Rcusp)<0.4
Rfold_s,Rcusp_s,phi0_s,phi1_s=Rfold_s[mask],Rcusp_s[mask],phi0_s[mask],phi1_s[mask]
#Rcusp_s=Rcusp[np.abs(Rcusp)<1.0]

#Rfold=Rfold[Rfold<0.5]
#Rfold_s=Rfold_s[Rfold_s<0.5]

print Rfold.size, Rfold_s.size,phi0.size,phi1_s.size
#### 

## ---- probability contour
## ellip
edge = np.linspace(40,120,9)
bmark = np.empty(Rfold.size)
bmark.fill(-1)

Rc_mid = np.empty(edge.size-1)
Rc_mid.fill(-1)

Ra_mid = np.empty(edge.size-1)
Ra_mid.fill(-1)
Rb_mid = np.empty(edge.size-1)
Rb_mid.fill(-1)


for j in range(edge.size-1):
	bmark[np.logical_and(phi0>=edge[j],phi0<edge[j+1])] = j

	Rbin = np.abs(Rcusp[bmark == j])
	Rc_mid[j] = np.median(Rbin)

	Ra_bin = Rbin[Rbin>Rc_mid[j]]
	Ra_mid[j] = np.median(Ra_bin)
	Rb_bin = Rbin[Rbin<Rc_mid[j]]
	Rb_mid[j] = np.median(Rb_bin)

edge2 = np.linspace(0,60,7)
Rf_mid = np.empty(edge2.size-1)
Rf_mid.fill(-1)
bmark.fill(-1)

for j in range(edge2.size-1):
	bmark[np.logical_and(phi1>=edge2[j],phi1<edge2[j+1])] = j

	Rf_mid[j] = np.median(np.abs(Rfold[bmark == j]))

## disk

#edge = np.linspace(40,120,9)
bmark = np.empty(Rfold_s.size)
bmark.fill(-1)

Rc_mid_s = np.empty(edge.size-1)
Rc_mid_s.fill(-1)


for j in range(edge.size-1):
	bmark[np.logical_and(phi0_s>=edge[j],phi0_s<edge[j+1])] = j

	#Rf_mid[j] = np.median(np.abs(Rfold[bmark == j]))
	Rc_mid_s[j] = np.median(Rcusp_s[bmark == j])

#edge2 = np.linspace(0,60,7)
Rf_mid_s = np.empty(edge2.size-1)
Rf_mid_s.fill(-1)
bmark.fill(-1)

for j in range(edge2.size-1):
	bmark[np.logical_and(phi1_s>=edge2[j],phi1_s<edge2[j+1])] = j

	Rf_mid_s[j] = np.median(np.abs(Rfold_s[bmark == j]))


#print bmark[:20]
#print phi0[:20].astype(int)
#print Rf_mid

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


#H,xbin,ybin = np.histogram2d(phi0_s,np.abs(Rcusp_s),weights=w_Rfold_s,bins=(np.linspace(0,180,40),np.linspace(0,0.4,40)))
H,xbin,ybin = np.histogram2d(phi1,np.abs(Rfold),weights=w_Rfold,bins=(np.linspace(0,60,40),np.linspace(0,0.4,40)))
H=H.T
fig1 = plt.figure()
plt.imshow(H,interpolation='nearest',origin='low',extent=[xbin[0],xbin[-1],ybin[0],ybin[-1]])
plt.gca().set_aspect(60/0.4)


#plt.scatter(edge[1:-1],Rc_mid[1:],color='b',label='elliptical 50%')
#plt.scatter(edge[1:-1],Ra_mid[1:],color='b',marker='*',label='elliptical 25%')
#plt.scatter(edge[1:-1],Rc_mid_s[1:],color='r',label='disk')

#plt.title('disks half area (pt~750)')
#plt.title('disks point source (pt~2000)')
#plt.title('elliptical half area (pt~750)')
plt.title('elliptical point source (pt~2000)')
#plt.ylabel('|Rcusp|')
plt.ylabel('|Rfold|')
#plt.ylabel('number counts')
#plt.xlabel('delta phi')
plt.xlabel('phi 1')
plt.legend()
plt.colorbar()
#plt.show()
plt.savefig('../../data/glamer/glamer_elp_pt_fold.png')

