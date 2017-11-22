import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.mplot3d import Axes3D
from astropy.cosmology import WMAP9 as cosmo

real_list_good = np.array([252,  493, 1436, 1916, 2027, 2108, 2179, 2458, 2742, 4808, 4814,
       6862, 7285, 7618, 8968, 9411])
#real_list_good = np.array([1945])
real_list_all = np.arange(10000)
mask = ~np.in1d(real_list_all,real_list_good)
real_list_select = real_list_all[mask]

real_list_bad = np.random.choice(real_list_select,len(real_list_good),replace=False)

real_list_good = real_list_good.astype(str)
real_list_bad = real_list_bad.astype(str)
#real_list = np.array([138]).astype(str)
print len(real_list_good)

re = 0.75
dist_cut_los = 0.3
dist_cut_sub = 0.12
area_los,area_sub = np.pi*dist_cut_los**2,np.pi*dist_cut_sub**2

zl,zs = 0.34,3.62
Ds, Dl = cosmo.angular_diameter_distance(zs),cosmo.angular_diameter_distance(zl)

avg_kappa_good = []
avg_kappa_bad = []

kappa0_good_sub,kappa1_good_sub,kappa2_good_sub = [],[],[]
kappa0_good_los,kappa1_good_los,kappa2_good_los = [],[],[]
kappa0_bad_sub,kappa1_bad_sub,kappa2_bad_sub = [],[],[]
kappa0_bad_los,kappa1_bad_los,kappa2_bad_los = [],[],[]

def cal_beta(z_los):

	beta = np.zeros(len(z_los))
	for i in range(len(z_los)):

		zi = np.array([zl,z_los[i]])
		z1,z2 = np.min(zi),np.max(zi)

		D1,D2 = cosmo.angular_diameter_distance(z1),cosmo.angular_diameter_distance(z2)
		

		beta[i] = (D2-D1)*Ds/D2/(Ds-D1)
	return beta

def dist_ratio(z_los):

	Dr = np.zeros(len(z_los))
	D_los = cosmo.angular_diameter_distance(z_los)

	for i in range(len(z_los)):

		if z_los[i]<zl:
			Dr[i] = D_los[i]/Dl

		else:
			Dr[i] = (Ds-D_los[i])/(Ds-Dl)

	return	Dr

img_x = np.array([0.0,0.38925,-0.33388,0.95065])
img_y = np.array([0.0,0.31998,-0.74771,-0.80215])


for  realid in real_list_good:
	los_table = np.loadtxt('/Volumes/sting_1/subs/los_00/los'+realid+'.txt')
	sub_table = np.loadtxt('/Volumes/sting_1/subs/real_01/real'+realid+'.txt')
	
	## ---- los
	z=los_table[:,0]
	mass = los_table[:,1]
	x = los_table[:,2]+7.65020767e-01
	y = los_table[:,3]-6.73284706e-01

	beta = np.abs(cal_beta(z))
	Dr = dist_ratio(z)

	#x,y = x/Dr,y/Dr

	dist0 = np.sqrt((x-img_x[0])**2+(y-img_y[0])**2)
	dist1 = np.sqrt((x-img_x[1])**2+(y-img_y[1])**2)
	dist2 = np.sqrt((x-img_x[2])**2+(y-img_y[2])**2)

	mask0 = dist0<dist_cut_los
	kappa0 = np.sum(mass[mask0]*beta[mask0])/area_los
	#print kappa0
	#kappa0 = np.sum(mass[mask0])/area

	mask1 = dist1<dist_cut_los
	kappa1 = np.sum(mass[mask1]*beta[mask1])/area_los
	#kappa1 = np.sum(mass[mask1])/area

	mask2 = dist2<dist_cut_los
	kappa2 = np.sum(mass[mask2]*beta[mask2])/area_los
	#kappa2 = np.sum(mass[mask2])/area

	kappa0_good_los.append(kappa0)
	kappa1_good_los.append(kappa1)
	kappa2_good_los.append(kappa2)


	## --- sub
	mass = sub_table[:,0]
	x = sub_table[:,1]+7.65020767e-01
	y = sub_table[:,2]-6.73284706e-01

	dist0 = np.sqrt((x-img_x[0])**2+(y-img_y[0])**2)
	dist1 = np.sqrt((x-img_x[1])**2+(y-img_y[1])**2)
	dist2 = np.sqrt((x-img_x[2])**2+(y-img_y[2])**2)

	mask0 = dist0<dist_cut_sub
	kappa0 = np.sum(mass[mask0])/area_sub

	mask1 = dist1<dist_cut_sub
	kappa1 = np.sum(mass[mask1])/area_sub

	mask2 = dist2<dist_cut_sub
	kappa2 = np.sum(mass[mask2])/area_sub

	#print np.log10(kappa0),np.log10(kappa1),np.log10(kappa2)
	avg_kappa_good.append(np.log10(kappa0+kappa1+kappa2))

	kappa0_good_sub.append(kappa0)
	kappa1_good_sub.append(kappa1)
	kappa2_good_sub.append(kappa2)

	kappa0_good_sub[kappa0_good_sub<0] = 0.0
	kappa1_good_sub[kappa0_good_sub<0] = 0.0
	kappa2_good_sub[kappa0_good_sub<0] = 0.0

	'''
	## --- surface density bin
	plt.hist2d(x,y,bins = 7, weights=mass,norm=LogNorm()) # 10
	
	#plt.scatter(x,y,marker='o',facecolor='none',edgecolor='b',s=mass/1e6)
	plt.colorbar()
	plt.clim(5e6,5e7)
	plt.scatter(img_x,img_y,marker='*',color='r',s=50)
	plt.title('chi2 < 10')
	plt.gca().set_aspect('equal')
	#plt.show()
	plt.savefig('/Volumes/sting_1/subs/real_01/kappa_real'+realid+'.png')
	plt.clf()
	'''

for  realid in real_list_bad:
	los_table = np.loadtxt('/Volumes/sting_1/subs/los_00/los'+realid+'.txt')
	sub_table = np.loadtxt('/Volumes/sting_1/subs/real_01/real'+realid+'.txt')
	
	## ---- los
	z=los_table[:,0]
	mass = los_table[:,1]
	x = los_table[:,2]+7.65020767e-01
	y = los_table[:,3]-6.73284706e-01

	beta = np.abs(cal_beta(z))
	Dr = dist_ratio(z)

	#x,y = x/Dr,y/Dr

	dist0 = np.sqrt((x-img_x[0])**2+(y-img_y[0])**2)
	dist1 = np.sqrt((x-img_x[1])**2+(y-img_y[1])**2)
	dist2 = np.sqrt((x-img_x[2])**2+(y-img_y[2])**2)

	mask0 = dist0<dist_cut_los
	kappa0 = np.sum(mass[mask0]*beta[mask0])/area_los
	#kappa0 = np.sum(mass[mask0])/area

	mask1 = dist1<dist_cut_los
	kappa1 = np.sum(mass[mask1]*beta[mask1])/area_los
	#kappa1 = np.sum(mass[mask1])/area

	mask2 = dist2<dist_cut_los
	kappa2 = np.sum(mass[mask2]*beta[mask2])/area_los
	#kappa2 = np.sum(mass[mask2])/area

	kappa0_bad_los.append(kappa0)
	kappa1_bad_los.append(kappa1)
	kappa2_bad_los.append(kappa2)

	## --- sub
	mass = sub_table[:,0]
	x = sub_table[:,1]+7.65020767e-01
	y = sub_table[:,2]-6.73284706e-01

	dist0 = np.sqrt((x-img_x[0])**2+(y-img_y[0])**2)
	dist1 = np.sqrt((x-img_x[1])**2+(y-img_y[1])**2)
	dist2 = np.sqrt((x-img_x[2])**2+(y-img_y[2])**2)

	mask0 = dist0<dist_cut_sub
	kappa0 = np.sum(mass[mask0])/area_sub

	mask1 = dist1<dist_cut_sub
	kappa1 = np.sum(mass[mask1])/area_sub

	mask2 = dist2<dist_cut_sub
	kappa2 = np.sum(mass[mask2])/area_sub

	#print np.log10(kappa0),np.log10(kappa1),np.log10(kappa2)
	avg_kappa_bad.append(np.log10(kappa0+kappa1+kappa2))

	kappa0_bad_sub.append(kappa0)
	kappa1_bad_sub.append(kappa1)
	kappa2_bad_sub.append(kappa2)

	kappa0_bad_sub[kappa0_bad_sub<0] = 0.0
	kappa1_bad_sub[kappa0_bad_sub<0] = 0.0
	kappa2_bad_sub[kappa0_bad_sub<0] = 0.0

#### -------

kappa0_bad_los,kappa1_bad_los,kappa2_bad_los = np.array(kappa0_bad_los),np.array(kappa1_bad_los),np.array(kappa2_bad_los)
kappa0_bad_sub,kappa1_bad_sub,kappa2_bad_sub = np.array(kappa0_bad_sub),np.array(kappa1_bad_sub),np.array(kappa2_bad_sub)
#kappa0_bad,kappa1_bad,kappa2_bad = kappa0_bad[kappa0_bad>0],kappa1_bad[kappa1_bad>0],kappa2_bad[kappa2_bad>0]

kappa0_good_los,kappa1_good_los,kappa2_good_los = np.array(kappa0_good_los),np.array(kappa1_good_los),np.array(kappa2_good_los)
kappa0_good_sub,kappa1_good_sub,kappa2_good_sub = np.array(kappa0_good_sub),np.array(kappa1_good_sub),np.array(kappa2_good_sub)
#kappa0_good,kappa1_good,kappa2_good = kappa0_good[kappa0_good>0],kappa1_good[kappa1_good>0],kappa2_good[kappa2_good>0]

print np.log10(kappa1_good_los)
print np.log10(kappa1_good_sub)

print np.log10(kappa1_bad_los)
print np.log10(kappa1_bad_sub)
#print np.log10(kappa1_bad_sub+kappa1_bad_los)

plt.hist(np.log10(kappa1_bad_sub+1.5*kappa1_bad_los),color='r',label='chi2 ~ smooth')
plt.hist(np.log10(kappa1_good_sub+1.5*kappa1_good_los),label='chi2 < 10',alpha=0.5)

#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#ax.scatter(np.log10(kappa0_good_sub+kappa0_good_los),np.log10(kappa1_good_sub+kappa1_good_los),np.log10(kappa2_good_sub+kappa2_good_los))
#ax.scatter(np.log10(kappa0_bad_sub+kappa0_bad_los),np.log10(kappa1_bad_sub+kappa1_bad_los),np.log10(kappa2_bad_sub+kappa2_bad_los),marker='^',c='r')

#plt.scatter(np.log10(kappa0_good_sub+kappa0_good_los),np.log10(kappa1_good_sub+kappa1_good_los))
#plt.scatter(np.log10(kappa0_bad_sub+kappa0_bad_los),np.log10(kappa1_bad_sub+kappa1_bad_los),c='r')
#plt.legend(loc=2)
#plt.xlim(7.5,9.2)
#plt.xlabel(' imgB surface mass density (log(Msun)/arcsec^2)')
#plt.ylabel(' imgA surface mass density (log(Msun)/arcsec^2)')
#plt.ylabel('number counts')
#plt.title('image A, sub+LOS')
#ax.set_xlabel('image B')
#ax.set_ylabel('image A')
#ax.set_zlabel('image C')

plt.show()

#plt.savefig('los00_imgC.png')