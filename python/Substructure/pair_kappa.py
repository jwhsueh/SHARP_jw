import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.mplot3d import Axes3D
from astropy.cosmology import WMAP9 as cosmo

pair_file = "/Volumes/sting_1/subs/B1422_0100_pair.txt"
idx_list = np.loadtxt(pair_file,dtype='int')
sub_list,los_list = idx_list[:,0],idx_list[:,1]
real_list_good = np.array([ 181,    820,   2269,   4630,   4636,   5857,   8840,   8964,
         9101,   9505,  10760,  11336,  11560,  11719,  12644,  13828,
        13925,  17008,  17622,  17780,  17865,  20532,  20901,  21066,
        22741,  22872,  23418,  23732,  23759,  23878,  24890,  25633,
        25751,  26937,  27685,  27897,  28663,  29819,  30098,  30434,
        30521,  30649,  32354,  33980,  35764,  36076,  36130,  36221,
        36536,  39244,  39496,  40421,  40619,  40980,  41951,  43676,
        43890,  44696,  45163,  45922,  46039,  47238,  47576,  49580,
        51452,  52735,  53011,  53050,  53688,  56944,  57281,  58300,
        58556,  59020,  59731,  60702,  61103,  61386,  64223,  65308,
        65508,  66320,  66589,  66673,  68026,  68851,  70109,  70253,
        70296,  70842,  71242,  71624,  72939,  72949,  74247,  74437,
        76362,  77243,  77385,  78845,  79414,  79798,  80334,  81474,
        83238,  83662,  86139,  86269,  86530,  87129,  87806,  87816,
        89419,  90091,  90959,  92219,  93405,  93546,  94412,  94442,
        94628,  95171,  96896,  97023,  97030,  97202,  98584,  99054,
       100523, 100541, 100922, 101373, 105429, 105927, 107668, 107992,
       108539, 108796, 109519, 109930])
#real_list_good = np.array([1945])
real_list_all = np.arange(110000)
mask = ~np.in1d(real_list_all,real_list_good)
real_list_select = real_list_all[mask]

real_list_bad = np.random.choice(real_list_select,len(real_list_good),replace=False)

#real_list_good = real_list_good.astype(str)
#real_list_bad = real_list_bad.astype(str)
#real_list = np.array([138]).astype(str)
print len(real_list_good)

re = 0.75
dist_cut_los = 0.5
dist_cut_sub = 0.13
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


for  realidx in real_list_good:

	subid,losid = sub_list[realidx],los_list[realidx]
	los_table = np.loadtxt('/Volumes/sting_1/subs/los_00/los'+str(losid)+'.txt')
	sub_table = np.loadtxt('/Volumes/sting_1/subs/real_01/real'+str(subid)+'.txt')
	
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

for  realidx in real_list_bad:
	subid,losid = sub_list[realidx],los_list[realidx]
	los_table = np.loadtxt('/Volumes/sting_1/subs/los_00/los'+str(losid)+'.txt')
	sub_table = np.loadtxt('/Volumes/sting_1/subs/real_01/real'+str(subid)+'.txt')
	
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

plt.hist(np.log10(kappa1_bad_sub+kappa1_bad_los),color='r',label='chi2 ~ smooth')
plt.hist(np.log10(kappa1_good_sub+kappa1_good_los),label='chi2 < 10',alpha=0.5)

#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#ax.scatter(np.log10(kappa0_good_sub+kappa0_good_los),np.log10(kappa1_good_sub+kappa1_good_los),np.log10(kappa2_good_sub+kappa2_good_los))
#ax.scatter(np.log10(kappa0_bad_sub+kappa0_bad_los),np.log10(kappa1_bad_sub+kappa1_bad_los),np.log10(kappa2_bad_sub+kappa2_bad_los),marker='^',c='r')

#plt.scatter(np.log10(kappa0_good_sub+1.5*kappa0_good_los),np.log10(kappa1_good_sub+1.5*kappa1_good_los))
#plt.scatter(np.log10(kappa0_bad_sub+1.5*kappa0_bad_los),np.log10(kappa1_bad_sub+1.5*kappa1_bad_los),c='r',marker='x')
plt.legend(loc=2)
plt.xlim(7.0,9.5)
plt.xlabel('surface mass density (log(Msun)/arcsec^2)')
#plt.ylabel(' imgA surface mass density (log(Msun)/arcsec^2)')
plt.ylabel('number counts')
plt.title('image A, sub+LOS hist')
#ax.set_xlabel('image B')
#ax.set_ylabel('image A')
#ax.set_zlabel('image C')

#plt.show()

plt.savefig('0100pair_imgA_all.png')