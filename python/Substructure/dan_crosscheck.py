import numpy as np
import matplotlib.pyplot as plt

basePath = '../../data/illustris_1'
ssNumber = '99'

catalog = basePath+'/Galaxy_Lens'+ssNumber+'.dat'
#GalaxyID = np.loadtxt(catalog,dtype = 'int',unpack=True, usecols=[0])
Mass = np.loadtxt(catalog,dtype = 'float',unpack=True, usecols=[1])
morph = np.loadtxt(catalog,dtype = 'int',unpack=True, usecols=[13])

Re = np.loadtxt(catalog,dtype = 'float',unpack=True, usecols=[7,8,9])
DMfrac = np.loadtxt(catalog,dtype = 'float',unpack=True, usecols=[10,11,12])

theta = np.loadtxt(catalog,dtype = 'float',unpack=True, usecols=[4,5,6])

# reduce to 1d [x...y...z...]
Re,DMfrac,theta = np.ravel(Re),np.ravel(DMfrac),np.ravel(theta) 

# expend array

Mass = np.array([Mass,Mass,Mass]).flatten()
print Mass

## Dandan's pick

# three projections

morph = morph.flatten()
#print morph
morph = np.array([morph,morph,morph]).flatten()
cri = morph==1

Re_dan = Re[cri]
DMfrac_dan = DMfrac[cri]
theta_dan = theta[cri]
Mass_dan = Mass[cri]

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

print Mass_dan_ed.size,Mass_dan_fa.size

## Kinematics pick

df = np.loadtxt(catalog,dtype = 'float',unpack=True, usecols=[2]) # disk str frac
bf = np.loadtxt(catalog,dtype = 'float',unpack=True, usecols=[3])

df = np.array([df,df,df]).flatten()
bf = np.array([bf,bf,bf]).flatten()

cri_k1 = df==1
cri_k2 = bf==1

Re_df,Re_bf = Re[cri_k1],Re[cri_k2]
Mass_df,Mass_bf = Mass[cri_k1],Mass[cri_k2]
theta_df,theta_bf = theta[cri_k1],theta[cri_k2]

# Edge-on/Face-on
edge_df,edge_bf = theta_df == 1,theta_bf == 1
face_df,face_bf = theta_df == 0,theta_bf == 0

Re_df_ed, Re_df_fa = Re_df[edge_df],Re_df[face_df]
Re_bf_ed, Re_bf_fa = Re_bf[edge_bf],Re_bf[face_bf]

Mass_df_ed, Mass_df_fa = Mass_df[edge_df],Mass_df[face_df]
Mass_bf_ed, Mass_bf_fa = Mass_bf[edge_bf],Mass_bf[face_bf]

print Mass_df_ed.size,Mass_df_fa.size
print Mass_bf_ed.size,Mass_bf_fa.size

# histogram
se = np.linspace(0.2,1.0,20)
dot = []
for i in range(se.size-1):
	dot.append((se[i]+se[i+1])/2.)

All = np.histogram(Re,bins = se)[0].astype(float)
Dan = np.histogram(Re_dan,bins = se)[0].astype(float)
Dan_ed = np.histogram(Re_dan_ed,bins = se)[0].astype(float)
Dan_fa = np.histogram(Re_dan_fa,bins = se)[0].astype(float)

df = np.histogram(Re_df,bins = se)[0].astype(float)
df_ed = np.histogram(Re_df_ed,bins = se)[0].astype(float)
df_fa = np.histogram(Re_df_fa,bins = se)[0].astype(float)

bf = np.histogram(Re_bf,bins = se)[0].astype(float)
bf_ed = np.histogram(Re_bf_ed,bins = se)[0].astype(float)
bf_fa = np.histogram(Re_bf_fa,bins = se)[0].astype(float)

plt.plot(dot,df/All,color='k',label = 'Disk str > 40%')
plt.plot(dot,df_fa/All,color='g',label = 'Disk str Face-on')
plt.plot(dot,df_ed/All,color='r',label = 'Disk str Edge-on')

plt.plot(dot,bf/All,'k--',label = 'Bulge str < 60%')
plt.plot(dot,bf_fa/All,'g--',label = 'Bulge str Face-on')
plt.plot(dot,bf_ed/All,'r--',label = 'Bulge str Edge-on')

#plt.scatter(np.log10(Mass_dan_fa),DMfrac_dan_fa,edgecolor='b',facecolors = 'none',marker='o',label='Dandan Face-on')
#plt.scatter(np.log10(Mass_dan_ed),DMfrac_dan_ed,edgecolor='r',facecolors = 'none',marker='^',label='Dandan Edge-on')


plt.xlabel('Einstein Radius')
plt.ylabel('Galaxy count frac')
plt.title('Snapshot99: 3 projections')
#plt.legend(scatterpoints=1,loc =2)
plt.legend()
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