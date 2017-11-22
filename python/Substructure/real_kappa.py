import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.mplot3d import Axes3D

real_list_good = np.array([  1195,  2298,  5305,  6174, 10171, 10309, 15108, 15720, 16643,
       17303, 18985, 20236, 22148, 23844, 24860, 28933, 32147, 33646,
       33934, 34303, 35860])
#real_list_good =np.array([1196])
real_list_all = np.arange(35700)
mask = ~np.in1d(real_list_all,real_list_good)
real_list_select = real_list_all[mask]

real_list_bad = np.random.choice(real_list_select,len(real_list_good),replace=False)

real_list_good = real_list_good.astype(str)
real_list_bad = real_list_bad.astype(str)
#real_list = np.array([138]).astype(str)
print len(real_list_good)

lens='MG0414'
re = 1.11
# B1422
#img_x = np.array([0.0,0.38925,-0.33388,0.95065])
#img_y = np.array([0.0,0.31998,-0.74771,-0.80215])

img_x = np.array([0.5876,0.7208,0.0,-1.3608])
img_y = np.array([-1.9341,-1.5298,0.0,-1.6348])
x0,y0 =  -5.569528e-01, -1.349765e+00

dist_cut = 0.1
area = np.pi*dist_cut**2

avg_kappa_good = []
avg_kappa_bad = []

kappa0_good,kappa1_good,kappa2_good = [],[],[]
kappa0_bad,kappa1_bad,kappa2_bad = [],[],[]

for  realid in real_list_good:
	los_table = np.loadtxt('/Volumes/sting_1/subs/'+lens+'/real_01/real'+realid+'.txt')
	#z=los_table[:,0]
	mass = los_table[:,0]
	x = los_table[:,1]+x0
	y = los_table[:,2]+y0



	dist0 = np.sqrt((x-img_x[0])**2+(y-img_y[0])**2)
	dist1 = np.sqrt((x-img_x[1])**2+(y-img_y[1])**2)
	dist2 = np.sqrt((x-img_x[2])**2+(y-img_y[2])**2)

	mask0 = dist0<dist_cut
	kappa0 = np.sum(mass[mask0])/area

	mask1 = dist1<dist_cut
	kappa1 = np.sum(mass[mask1])/area

	mask2 = dist2<dist_cut
	kappa2 = np.sum(mass[mask2])/area

	#print np.log10(kappa0),np.log10(kappa1),np.log10(kappa2)
	avg_kappa_good.append(np.log10(kappa0+kappa1+kappa2))

	kappa0_good.append(np.log10(kappa0))
	kappa1_good.append(np.log10(kappa1))
	kappa2_good.append(np.log10(kappa2))

	
	## --- surface density bin
	plt.hist2d(x,y,bins = 7, weights=mass,norm=LogNorm()) # 10
	
	#plt.scatter(x,y,marker='o',facecolor='none',edgecolor='b',s=mass/1e6)
	plt.colorbar()
	plt.clim(5e6,5e7)
	plt.scatter(img_x,img_y,marker='*',color='r',s=50)
	plt.title('chi2 < 1.5')
	plt.gca().set_aspect('equal')
	#plt.show()
	#plt.savefig('/Volumes/sting_1/subs/'+lens+'/real_01/kappa_real'+realid+'.png')
	plt.savefig('kappa_real'+realid+'.png')
	plt.clf()
	
'''
for  realid in real_list_bad:
	los_table = np.loadtxt('/Volumes/sting_1/subs/'+lens+'/real_01/real'+realid+'.txt')
	#z=los_table[:,0]
	mass = los_table[:,0]
	x = los_table[:,1]+x0
	y = los_table[:,2]+y0

	dist0 = np.sqrt((x-img_x[0])**2+(y-img_y[0])**2)
	dist1 = np.sqrt((x-img_x[1])**2+(y-img_y[1])**2)
	dist2 = np.sqrt((x-img_x[2])**2+(y-img_y[2])**2)

	mask0 = dist0<dist_cut
	kappa0 = np.sum(mass[mask0])/area

	mask1 = dist1<dist_cut
	kappa1 = np.sum(mass[mask1])/area

	mask2 = dist2<dist_cut
	kappa2 = np.sum(mass[mask2])/area

	#print np.log10(kappa0),np.log10(kappa1),np.log10(kappa2)
	avg_kappa_bad.append(np.log10(kappa0+kappa1+kappa2))

	kappa0_bad.append(np.log10(kappa0))
	kappa1_bad.append(np.log10(kappa1))
	kappa2_bad.append(np.log10(kappa2))

kappa0_bad,kappa1_bad,kappa2_bad = np.array(kappa0_bad),np.array(kappa1_bad),np.array(kappa2_bad)
kappa0_bad,kappa1_bad,kappa2_bad = kappa0_bad[kappa0_bad>0],kappa1_bad[kappa1_bad>0],kappa2_bad[kappa2_bad>0]

kappa0_good,kappa1_good,kappa2_good = np.array(kappa0_good),np.array(kappa1_good),np.array(kappa2_good)
kappa0_good,kappa1_good,kappa2_good = kappa0_good[kappa0_good>0],kappa1_good[kappa1_good>0],kappa2_good[kappa2_good>0]
#print np.array(kappa1_bad)>0

plt.hist(kappa1_bad,color='r',label='chi2 ~ smooth')
plt.hist(kappa1_good,label='chi2 < 2',alpha=0.5)

#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#ax.scatter(kappa0_good,kappa1_good,kappa2_good)
#ax.scatter(kappa0_bad,kappa1_bad,kappa2_bad,marker='^',c='r')
plt.legend(loc=2)
#plt.xlim(7.5,9.2)
plt.xlabel('surface mass density (log(Msun)/arcsec^2)')
plt.ylabel('number counts')
#plt.title('image A1, 0.13" radii')
#ax.set_xlabel('image B')
#ax.set_ylabel('image A')
#ax.set_zlabel('image C')
plt.show()
'''
#plt.savefig('sub01_imgC.png')