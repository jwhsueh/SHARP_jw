import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.mplot3d import Axes3D

real_list_good = np.array([138,  589,  915, 1460, 3703, 3841, 3864, 4882,5106, 5461, 5768, 5928, 6077, 6862, 8023, 8850,
	10009, 10482, 10804, 11995, 13372, 13531, 13705, 14131, 16683,
       17805, 17883, 18353, 18552, 19045, 19432, 19487, 20071, 20904,
       21840, 21855, 22091, 22462, 23220, 23980, 24298, 24667, 24836,
       25236, 26243, 26498, 27232, 27887, 29179, 29536, 29725, 31358,
       31799, 32593, 33514, 33563, 36282, 37343, 37974, 38174, 38669,
       39395, 39774, 40651, 41386, 42610, 44993, 46338, 46448, 46633,
       48810, 49448, 49660, 50004, 50617, 50635, 50824, 51678, 52079,
       52270, 52453, 54197, 54975, 55189, 55679, 56100, 57181, 57226,
       58218, 58435,
       60307,  60579,  60865,  62020,  62206,  62240,  62270,  63813,
        64029,  64906,  65490,  65818,  66671,  67048,  67882,  68189,
        68640,  69485,  69746,  71152,  71405,  71985,  72505,  72601,
        72734,  72872,  72895,  73141,  74408,  76250,  77811,  78922,
        79451,  79547,  80862,  81078,  81538,  82263,  82783,  83122,
        83155,  84405,  84707,  85220,  86627,  87102,  87487,  88920,
        89148,  89836,  90006,  90478,  91199,  91784,  92589,  93326,
        94329,  94665,  95428,  95864,  96642,  96734,  97893,  98307,
        98824,  99533,  99710, 100783, 101003, 101545, 101860, 103929,
       105625, 106168, 106370, 106739, 106837, 106932, 108794, 109111,
       109677, 109816])
real_list_all = np.arange(110000)
mask = ~np.in1d(real_list_all,real_list_good)
real_list_select = real_list_all[mask]

real_list_bad = np.random.choice(real_list_select,len(real_list_good),replace=False)

real_list_good = real_list_good.astype(str)
real_list_bad = real_list_bad.astype(str)
#real_list = np.array([138]).astype(str)
print len(real_list_good)

re = 0.75
dist_cut = 0.12
area = np.pi*dist_cut**2

avg_kappa_good = []
avg_kappa_bad = []

kappa0_good,kappa1_good,kappa2_good = [],[],[]
kappa0_bad,kappa1_bad,kappa2_bad = [],[],[]

for  realid in real_list_good:
	los_table = np.loadtxt('/Volumes/sting_1/subs/real_01/real'+realid+'.txt')
	#z=los_table[:,0]
	mass = los_table[:,0]
	x = los_table[:,1]+7.65020767e-01
	y = los_table[:,2]-6.73284706e-01

	img_x = np.array([0.0,0.38925,-0.33388,0.95065])
	img_y = np.array([0.0,0.31998,-0.74771,-0.80215])

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
	los_table = np.loadtxt('/Volumes/sting_1/subs/real_01/real'+realid+'.txt')
	#z=los_table[:,0]
	mass = los_table[:,0]
	x = los_table[:,1]+7.65020767e-01
	y = los_table[:,2]-6.73284706e-01

	img_x = np.array([0.0,0.38925,-0.33388,0.95065])
	img_y = np.array([0.0,0.31998,-0.74771,-0.80215])

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

plt.hist(kappa2_bad,color='r',label='chi2 ~ smooth')
plt.hist(kappa2_good,label='chi2 < 10',alpha=0.5)

#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#ax.scatter(kappa0_good,kappa1_good,kappa2_good)
#ax.scatter(kappa0_bad,kappa1_bad,kappa2_bad,marker='^',c='r')
plt.legend(loc=2)
#plt.xlim(7.5,9.2)
plt.xlabel('surface mass density (log(Msun)/arcsec^2)')
plt.ylabel('number counts')
plt.title('image C, 0.12" radii')
#ax.set_xlabel('image B')
#ax.set_ylabel('image A')
#ax.set_zlabel('image C')
#plt.show()

plt.savefig('sub01_imgC.png')