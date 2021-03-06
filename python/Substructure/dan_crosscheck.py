import numpy as np
import matplotlib.pyplot as plt
import DistanceTool as distance

basePath = '../../data/illustris_1'
ssNumber = '99'
z = 0.6

## change
catalog = basePath+'/Galaxy_Lens'+ssNumber+'_sig.dat'
#GalaxyID = np.loadtxt(catalog,dtype = 'int',unpack=True, usecols=[0])
Mass = np.loadtxt(catalog,dtype = 'float',unpack=True, usecols=[1])
str_ms = np.loadtxt(catalog,dtype = 'float',unpack=True, usecols=[2])
sigma = np.loadtxt(catalog,dtype = 'float',unpack=True, usecols=[3])
morph = np.loadtxt(catalog,dtype = 'int',unpack=True, usecols=[18])

Re = np.loadtxt(catalog,dtype = 'float',unpack=True, usecols=[9,10,11])
DMfrac = np.loadtxt(catalog,dtype = 'float',unpack=True, usecols=[12,13,14])

theta = np.loadtxt(catalog,dtype = 'float',unpack=True, usecols=[6,7,8])

relax = np.loadtxt(catalog,dtype = 'int',unpack=True, usecols=[19])
V_mag = np.loadtxt(catalog,dtype = 'float',unpack=True, usecols=[20])

SB = np.loadtxt(catalog,dtype = 'float',unpack=True, usecols=[15,16,17])

# reduce to 1d [x...y...z...]
Re,DMfrac,theta = np.ravel(Re),np.ravel(DMfrac),np.ravel(theta)
SB = np.ravel(SB) 

# expend array

Mass = np.array([Mass,Mass,Mass]).flatten()
str_ms = np.array([str_ms,str_ms,str_ms]).flatten()
sigma = np.array([sigma,sigma,sigma]).flatten()
relax = np.array([relax,relax,relax]).flatten()
V_mag = np.array([V_mag,V_mag,V_mag]).flatten()
#print Mass
mag_cut = 24.
SB_cut = 24


## magnitude cut for all
#cri0 = V_mag<mag_cut
cri0 = SB<SB_cut
Re_mag = Re[cri0]

## Dandan's pick

# three projections

morph = morph.flatten()
#print morph
morph = np.array([morph,morph,morph]).flatten()
cri = morph==1

Dan_relax = relax[cri]
Re_dan = Re[cri]
DMfrac_dan = DMfrac[cri]
theta_dan = theta[cri]
Mass_dan = Mass[cri]
str_dan = str_ms[cri]
sig_dan = sigma[cri]
Vmag_dan = V_mag[cri]
SB_dan = SB[cri]

## relax
cri_relax = Dan_relax ==0
Re_dan = Re_dan[cri_relax]
Vmag_dan = Vmag_dan[cri_relax]
SB_dan = SB_dan[cri_relax]

# mag cut
#cri_V = Vmag_dan<mag_cut
cri_V = SB_dan <SB_cut
Re_dan = Re_dan[cri_V]


'''
## surface brightness cut
cri2 = brightness_dan<25.5
#print cri2

Re_dan = Re_dan[cri2]
#print len(Re_dan)
'''
'''
# Edge-on/Face-on
edge = theta_dan == 1
face = theta_dan == 0

Re_dan_ed = Re_dan[edge]
Re_dan_fa = Re_dan[face]

DMfrac_dan_ed = DMfrac_dan[edge]
DMfrac_dan_fa = DMfrac_dan[face]

Mass_dan_ed = Mass_dan[edge]
#print Mass_dan_ed
Mass_dan_fa = Mass_dan[face]

str_dan_ed = str_dan[edge]
str_dan_fa = str_dan[face]

sig_dan_ed = sig_dan[edge]
sig_dan_fa = sig_dan[face]

print Mass_dan_ed.size,Mass_dan_fa.size
'''
## Kinematics pick

df = np.loadtxt(catalog,dtype = 'float',unpack=True, usecols=[4]) # disk str frac
bf = np.loadtxt(catalog,dtype = 'float',unpack=True, usecols=[5])

df = np.array([df,df,df]).flatten()
bf = np.array([bf,bf,bf]).flatten()

cri_k1 = df==1.
cri_k2 = bf==1.

Re_df,Re_bf = Re[cri_k1],Re[cri_k2]
#Re_bf =Re[cri_k2]
#Mass_df,Mass_bf = Mass[cri_k1],Mass[cri_k2]
Mass_bf = Mass[cri_k2]
#theta_df,theta_bf = theta[cri_k1],theta[cri_k2]
theta_bf = theta[cri_k2]

relax_df,relax_bf = relax[cri_k1],relax[cri_k2]
Vmag_df,Vmag_bf = V_mag[cri_k1],V_mag[cri_k2]
SB_df,SB_bf = SB[cri_k1],SB[cri_k2]

## relax

crik1_relax,crik2_relax = relax_df==0,relax_bf==0
Re_df,Re_bf = Re_df[crik1_relax],Re_bf[crik2_relax]
Vmag_df,Vmag_bf = Vmag_df[crik1_relax],Vmag_bf[crik2_relax]
SB_df,SB_bf = SB_df[crik1_relax],SB_bf[crik2_relax]

## mag cut

#crik1_mag,crik2_mag = Vmag_df<mag_cut,Vmag_bf<mag_cut
crik1_mag,crik2_mag = SB_df<SB_cut,SB_bf<SB_cut
Re_df,Re_bf = Re_df[crik1_mag],Re_bf[crik2_mag]
# add phot cut

#brightness_df,brightness_bf = brightness[cri_k1],brightness[cri_k2]
#cri_k11,cri_k22 = brightness_df < 23.5, brightness_bf < 23.5

#Re_df,Re_bf = Re_df[cri_k11],Re_bf[cri_k22]

'''
# Edge-on/Face-on
#edge_df,edge_bf = theta_df == 1,theta_bf == 1
#face_df,face_bf = theta_df == 0,theta_bf == 0

edge_bf = theta_bf == 1
face_bf = theta_bf == 0

#Re_df_ed, Re_df_fa = Re_df[edge_df],Re_df[face_df]
Re_bf_ed, Re_bf_fa = Re_bf[edge_bf],Re_bf[face_bf]

#Mass_df_ed, Mass_df_fa = Mass_df[edge_df],Mass_df[face_df]
#Mass_bf_ed, Mass_bf_fa = Mass_bf[edge_bf],Mass_bf[face_bf]

#print Mass_df_ed.size,Mass_df_fa.size
#print Mass_bf_ed.size,Mass_bf_fa.size
'''

## ---- additional surface brightness cut
'''
for i in range(Re_dan.size):
	if SB_dan[i]>SB_cut:
		list(Re_dan).pop(i)
'''
for i in range(Re_df.size):
	if SB_df[i]>SB_cut:
		list(Re_df).pop(i)

for i in range(Re_bf.size):
	if SB_bf[i]>SB_cut:
		list(Re_bf).pop(i)


# histogram
se = np.linspace(0.0,1.5,15)
dot = []
for i in range(se.size-1):
	dot.append((se[i]+se[i+1])/2.)

All = np.histogram(Re_mag,bins = se)[0].astype(float)
#All = np.histogram(Re,bins = se)[0].astype(float)
Dan = np.histogram(Re_dan,bins = se)[0].astype(float)
#Dan_ed = np.histogram(Re_dan_ed,bins = se)[0].astype(float)
#Dan_fa = np.histogram(Re_dan_fa,bins = se)[0].astype(float)

df = np.histogram(Re_df,bins = se)[0].astype(float)
#df_ed = np.histogram(Re_df_ed,bins = se)[0].astype(float)
#df_fa = np.histogram(Re_df_fa,bins = se)[0].astype(float)

bf = np.histogram(Re_bf,bins = se)[0].astype(float)
#bf_ed = np.histogram(Re_bf_ed,bins = se)[0].astype(float)
#bf_fa = np.histogram(Re_bf_fa,bins = se)[0].astype(float)



plt.plot(dot,Dan/All,color='k',label = 'morphology',linestyle='steps')
#plt.bar(np.array(dot)-0.03,Dan/All,edgecolor='k',facecolor='g',width=(se[1]-se[0])/5.,label = 'Morphology')

#plt.plot(dot,Dan_fa/All,color='g',label = 'morphology Face-on')
#plt.plot(dot,Dan_ed/All,color='r',label = 'morphology Edge-on')

plt.plot(dot,bf/All,'b-.',label = 'Bulge Star Fraction < 60%',linestyle = 'steps')
#plt.bar(dot,bf/All,edgecolor='k',facecolor='b',width=(se[1]-se[0])/5.,label = 'Bulge Star Fraction < 60%')

#plt.plot(dot,bf_fa/All,'g--',label = 'Bulge str Face-on')
#plt.plot(dot,bf_ed/All,'r--',label = 'Bulge str Edge-on')

plt.plot(dot,df/All,'r-.',label = 'Disk Star Fraction > 40%',linestyle = 'steps')
#plt.bar(np.array(dot)+0.03,df/All,edgecolor='k',facecolor='r',width=(se[1]-se[0])/5.,label = 'Disk Star Fraction > 40%')


#plt.scatter(Re_bf,SB_bf,edgecolor='r',facecolors = 'r',marker='^',label='Bulge Star Fraction < 60%')
#plt.scatter(Re_dan,SB_dan,edgecolor='b',facecolors = 'none',marker='o',label='Morphology')
#plt.scatter(Re_df,SB_df,edgecolor='k',facecolors = 'k',marker='v',label='Disk Star Fraction > 40%')
#plt.scatter(Re,V_mag,edgecolor='k',facecolors = 'none',marker='o',label='All galaxy')


plt.xlabel('Einstein Radius (arc sec)')
#plt.ylabel('Disk System Fraction')
plt.ylabel('r-band surface brightness')
plt.title('Effective surface brightness of ExpDisk')
plt.legend(scatterpoints=1,loc =1)
#plt.legend(loc = 1)
plt.xlim(0.2,1.2)
#plt.ylim(18,32)
plt.show()

'''
plt.hist(Re_dan_fa,bins,color='g',label = 'Dandan face-on',histtype='step')
plt.hist(Re_dan,bins,color='k',label = 'Dandan all',histtype='step')
plt.hist(Re_dan_ed,bins,color = 'r',label = 'Dandan Edge-on',histtype='step')


plt.hist(Galaxy_DMfrac,bins,label = 'All galaxy',histtype='step')
plt.hist(disk_DMfrac,bins,color = 'g',label = 'bulge Star frac<0.6',histtype='step')
plt.hist(DanDMfrac_mc,bins,color = 'r', label='Dandan morphology pick',histtype='step')

plt.hist(disk_DMfrac2,bins,color = 'k',label = 'disk Star frac>0.4',histtype='step')
plt.legend(loc =1)
plt.xlabel('DM frac w/i Re')
plt.ylabel('galaxy count')
plt.title('Snapshot 99')
plt.show()
'''

'''
plt.hist(disk_ms,bins,label = 'Bulge Star frac<0.6',)
plt.hist(Dan_ms,bins,color = 'r', label='Dandan morphology pick')
plt.hist(disk_ms2,bins,color = 'g',label = 'Disk Star frac>0.4',)
plt.legend()
plt.xlabel('M_sun')
plt.ylabel('galaxy count')
plt.title('Snapshot 99')
plt.show()
'''

'''

plt.ylim(30,16)
plt.plot(np.log10(Dan_ms),Dan_V,'+',label = 'morphology (Dandan)')
plt.plot(np.log10(disk_ms),disk_V,'r+',label = 'bulge star fraction<0.6')

#plt.plot(o_mass,o_frac,'o',mec ='g',mfc = 'none',label = 'overlap (beta<0.2)')
plt.xlabel('log10(Mass)')
plt.ylabel('V-band mag')
plt.legend()

plt.show()
'''