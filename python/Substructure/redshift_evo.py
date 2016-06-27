import numpy as np
import matplotlib.pyplot as plt

basePath = '../../data/illustris_1'

ssNumber = ['85','99','120']

All = []
disk_frac = []
bulge_frac = []
dot_z = []
for i in range(len(ssNumber)):
	catalog = basePath+'/Galaxy_Lens'+ssNumber[i]+'_sig.dat'

	Re_proj = np.loadtxt(catalog,dtype = 'float',unpack=True, usecols=[12,13,14])
	df = np.loadtxt(catalog,dtype = 'float',unpack=True, usecols=[4])
	bf = np.loadtxt(catalog,dtype = 'float',unpack=True, usecols=[5])

	# reduce to 1d [x...y...z...]

	Re = np.ravel(Re_proj)

	df = np.array([df,df,df]).flatten()
	bf = np.array([bf,bf,bf]).flatten()
	#disk_frac.append(df)
	#bulge_frac.append(bf)

	## plotting

	# criteria flag

	cri_k1 = df==1
	cri_k2 = bf==1

	Re_df,Re_bf = Re[cri_k1],Re[cri_k2]

	# histogram
	se = np.linspace(0.2,1.2,20)
	dot = []
	for i in range(se.size-1):
		dot.append((se[i]+se[i+1])/2.)

	dot_z.append(dot)

	All.append(np.histogram(Re,bins = se)[0].astype(float))
	disk_frac.append(np.histogram(Re_df,bins = se)[0].astype(float))
	bulge_frac.append(np.histogram(Re_bf,bins = se)[0].astype(float))


## plot

plt.plot(dot_z[0],disk_frac[0]/All[0],color='k',label = 'z = 1.0')
plt.plot(dot_z[1],disk_frac[1]/All[1],color='b',label = 'z = 0.6')
plt.plot(dot_z[2],disk_frac[2]/All[2],color='r',label = 'z = 0.2')

plt.legend()

plt.xlabel('Einstein Radius')
plt.ylabel('Galaxy Fraction')
plt.title('Redshift comparison: Disk star frac > 40%')

plt.show()



