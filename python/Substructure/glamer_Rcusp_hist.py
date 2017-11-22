import numpy as np
import matplotlib.pyplot as plt

path='/Volumes/sting_1/data/'
sie_path='/Users/jwhsueh/Documents/SHARP_jw/data/sub_gravlens/snap99_elp/'
sn_path = '/Volumes/sting_1/snap99_result/snap99_111/' ## particle SIE
sn3_path = '/Volumes/sting_1/snap99_result/snap99_222/'

list_file_sub = path+'snap99_sn2_Rcusp.txt'
list_file_sub2 = path+'snap99_sn3_Rcusp.txt'
list_file_e = path+'snap99_elp2_Rcusp.txt'
list_file_d = path+'snap99_tri2_Rcusp.txt'  # disc
#list_file_in = path+'snap99_edg_size.txt'
#size_tab = np.loadtxt(list_file_in)

file_list_sub = np.genfromtxt(list_file_sub,dtype='str')
file_list_sub2 = np.genfromtxt(list_file_sub2,dtype='str')
file_list_e = np.genfromtxt(list_file_e,dtype='str')
file_list_d = np.genfromtxt(list_file_d,dtype='str')
#theta = np.loadtxt(list_file_in)[:,2]
#file_list_sub = ['227953sie_p1']

Rfold_s,Rcusp_s,phi0_s,phi1_s  = np.empty(1),np.empty(1),np.empty(1),np.empty(1)
Rfold_e,Rcusp_e,phi0_e,phi1_e  = np.empty(1),np.empty(1),np.empty(1),np.empty(1)
Rfold_d,Rcusp_d,phi0_d,phi1_d  = np.empty(1),np.empty(1),np.empty(1),np.empty(1)

for file_one in file_list_sub:
	
	#table = np.loadtxt(path+'Rcusp_f/'+file_one+'_64_Rcusp_ga.txt')
	#table = np.loadtxt(sie_path+file_one+'_rcusp_cs.txt') # sie
	table = np.loadtxt(sn_path+file_one+'_Rcusp_ga.txt')
	#table = np.loadtxt(sie_path+file_one+'_rcusp.txt') # sie sn
	

	#mask = size_tab[:,3]>0.25
	#mask = mask[:(table[:,0]).size]
	#Rfold_s,Rcusp_s,phi0_s,phi1_s = np.append(Rfold_s,table[:,0][mask]),np.append(Rcusp_s,table[:,1][mask]),np.append(phi0_s,table[:,2][mask]),np.append(phi1_s,table[:,3][mask])
	Rfold_s,Rcusp_s,phi0_s,phi1_s = np.append(Rfold_s,table[:,0]),np.append(Rcusp_s,table[:,1]),np.append(phi0_s,table[:,2]),np.append(phi1_s,table[:,3])

	## -- add table 2 here

for file_one in file_list_sub2:
#	print file_one
	table = np.loadtxt(sn3_path+file_one+'_Rcusp_ga.txt')
	Rfold_s,Rcusp_s,phi0_s,phi1_s = np.append(Rfold_s,table[:,0]),np.append(Rcusp_s,table[:,1]),np.append(phi0_s,table[:,2]),np.append(phi1_s,table[:,3])

for file_one in file_list_e:
	table = np.loadtxt(path+'Rcusp_f/'+file_one+'_64_Rcusp_ga.txt')
	Rfold_e,Rcusp_e,phi0_e,phi1_e = np.append(Rfold_e,table[:,0]),np.append(Rcusp_e,table[:,1]),np.append(phi0_e,table[:,2]),np.append(phi1_e,table[:,3])

for file_one in file_list_d:
	table = np.loadtxt(path+'Rcusp_f/'+file_one+'_64_Rcusp_ga.txt')
	Rfold_d,Rcusp_d,phi0_d,phi1_d = np.append(Rfold_d,table[:,0]),np.append(Rcusp_d,table[:,1]),np.append(phi0_d,table[:,2]),np.append(phi1_d,table[:,3])


mask = np.abs(Rfold_s)<0.45
Rfold_s,Rcusp_s,phi0_s,phi1_s=Rfold_s[mask],Rcusp_s[mask],phi0_s[mask],phi1_s[mask]
mask = np.abs(Rcusp_s)<0.45
Rfold_s,Rcusp_s,phi0_s,phi1_s=Rfold_s[mask],Rcusp_s[mask],phi0_s[mask],phi1_s[mask]

mask = np.abs(Rfold_e)<0.7
Rfold_e,Rcusp_e,phi0_e,phi1_e=Rfold_e[mask],Rcusp_e[mask],phi0_e[mask],phi1_e[mask]
mask = np.abs(Rcusp_e)<0.7
Rfold_e,Rcusp_e,phi0_e,phi1_e=Rfold_e[mask],Rcusp_e[mask],phi0_e[mask],phi1_e[mask]

mask = np.abs(Rfold_d)<0.7
Rfold_d,Rcusp_d,phi0_d,phi1_d=Rfold_d[mask],Rcusp_d[mask],phi0_d[mask],phi1_d[mask]
mask = np.abs(Rcusp_d)<0.7
Rfold_d,Rcusp_d,phi0_d,phi1_d=Rfold_d[mask],Rcusp_d[mask],phi0_d[mask],phi1_d[mask]

#phi0_s = phi0_s-30
Rfold_s[phi1_s<30] = Rfold_s[phi1_s<30]/2

####
## histogram
mask = np.logical_and(phi1_s>5.,phi1_s<40)
Rfold_s = Rfold_s[mask]
#mask = np.logical_and(phi0_s>40.,phi0_s<110)
#Rcusp_s = Rcusp_s[mask]

#w = np.zeros(len(Rfold_s))
#w.fill(1./len(Rfold_s))
w = np.zeros(len(Rcusp_s))
w.fill(1./len(Rcusp_s))

#mask = np.logical_and(phi1_e>5.,phi1_e<40)
#Rfold_e = Rfold_e[mask]
mask = np.logical_and(phi0_e>4.,phi0_e<110.)
Rcusp_e = Rcusp_e[mask]

#w2 = np.zeros(len(Rfold_e))
#w2.fill(1./len(Rfold_e))
w2 = np.zeros(len(Rcusp_e))
w2.fill(1./len(Rcusp_e))

#mask = np.logical_and(phi1_d>5.,phi1_d<40)
#Rfold_d = Rfold_d[mask]
mask = np.logical_and(phi0_d>4.,phi0_d<110.)
Rcusp_d = Rcusp_d[mask]

#w3 = np.zeros(len(Rfold_d))
#w3.fill(1./len(Rfold_d))
w3 = np.zeros(len(Rcusp_d))
w3.fill(1./len(Rcusp_d))

b = np.arange(0.,0.7,0.05)

plt.hist(np.abs(Rcusp_s),weights=w,color='k',bins=b,label='Smooth',log='True',histtype='step',lw=2,ls='-.')
plt.hist(np.abs(Rcusp_e),weights=w2,color='r',bins=b,label='Elliptical',log='True',histtype='step',lw=2,ls='dashed')
plt.hist(np.abs(Rcusp_d),weights=w3,color='b',bins=b,label='Disc',log='True',histtype='step',lw=2)
#plt.hist(np.abs(Rfold_s),weights=w,color='k',alpha=0.5,bins=b,label='SIE',log='True')
#plt.hist(np.abs(Rfold_e),weights=w2,color='r',alpha=0.3,bins=b,label='ell',log='True')
plt.legend()
plt.xlabel('$|R_{cusp}|$')
plt.ylabel('Number fraction')
plt.ylim(0,0.6)
#plt.gca().set_aspect(0.5/0.7)
plt.show()
#plt.savefig('../../data/glamer/Rfold_hist_com_all.png')