import numpy as np
import matplotlib.pyplot as plt

caus=np.loadtxt('/Volumes/sting_1/subs/caustic_sie_0.txt')
#caus=np.loadtxt('/Users/jwhsueh/Documents/SHARP_jw/data/sub_gravlens/B1422_0.crit')
caus = np.degrees(caus)*3600
caus[:,0] = caus[:,0]
caus[:,1] = caus[:,1]

sx=3.906056671923252566e-01 
sy=-4.160411050175467640e-01

#mask2=crit2[:,2]==0
#mask=crit[:,2]==0
#caus_intv=np.sqrt((caus[1:,0]-caus[:-1,0])**2+(caus[1:,1]-caus[:-1,1])**2)
#print caus[1:,0]-caus[:-1,0]
#print min(caus_intv)

## converting


plt.figure(1)
#plt.subplot(211)7.453065220167933713e-01 7.676655103095577637e-01
#plt.scatter(caus[:,0],caus[:,1],color='r',marker='.',s=1)
plt.scatter((sx-7.676655103095577637e-01), (sy+6.754106455630485284e-01))
#plt.scatter(src_pt[10,0],src_pt[10,1],marker='*',s=50,color='b')
#plt.scatter(src2_pt[:,0],src_pt[:,1],marker='.',s=1)
plt.plot(caus[:,0],caus[:,1],color='r',label='no sub')
#plt.plot(caus2[:,0],caus2[:,1],color='k',label='w/ sub')

#plt.scatter(0.000225848, 0.000256699,color='r',marker='x')
#plt.scatter( 0.000225807, 0.000256587,color='k',marker='x')

#plt.scatter(0.000280112, -0.000279768,color='r',marker='o')
#plt.scatter(0.000280089, -0.000279743,color='k',marker='o')

#plt.scatter( 0.000279946+(0.000280089-0.000280112), -0.00027959+(-0.000279743+0.000279768),color='k',marker='x')
#plt.plot(crit[:,0][mask],crit[:,1][mask],color='k',label='w sub')
#plt.plot(crit[:,0],crit[:,1],color='k',label='no sub')

#plt.axis('equal')
#plt.ylim(-6.0*1e-6,-5.0*1e-6)
#plt.xlim(-4.0*1e-6,-2.0*1e-6)

#plt.ylim(-1.0*1e-5,4.0*1e-6)
#plt.xlim(-4.0*1e-6,1e-6)
#plt.title("NN=64 w sub (seperate)")
'''
plt.subplot(212)
#plt.plot(crit2[:,0][mask2],crit2[:,1][mask2],color='k',label='no sub')
plt.axis('equal')
plt.ylim(-6.0*1e-6,-5.0*1e-6)
plt.xlim(-4.0*1e-6,-2.0*1e-6)
plt.title("NN=64 w/o sub")
'''
plt.xlim(1.0,-1.0)
#plt.ylim(-2.0,1.0)
plt.gca().set_aspect('equal')
plt.show()

#plt.savefig('../../data/glamer/displace2_caustic.png')