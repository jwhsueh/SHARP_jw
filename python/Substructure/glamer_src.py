import numpy as np
import matplotlib.pyplot as plt

caus=np.loadtxt('/Users/jwhsueh/Documents/glamer/examples/ParticleExample/build/caustic_218142_p1_.txt')
caus2=np.loadtxt('/Users/jwhsueh/Documents/glamer/examples/ParticleExample/build/caustic_218142_p1_sh.txt')
#crit=np.loadtxt('/Users/jwhsueh/Documents/glamer/examples/ParticleExample/build/critical_proj3_281185sub_64.txt')
#crit2=np.loadtxt('/Users/jwhsueh/Documents/glamer/examples/ParticleExample/build/critical_proj3_281185_64.txt')
#src_pt=np.loadtxt('/Users/jwhsueh/Documents/glamer/examples/ParticleExample/build/src_position2.txt')

#mask2=crit2[:,2]==0
#mask=crit[:,2]==0
#caus_intv=np.sqrt((caus[1:,0]-caus[:-1,0])**2+(caus[1:,1]-caus[:-1,1])**2)
#print caus[1:,0]-caus[:-1,0]
#print min(caus_intv)

## converting


plt.figure(1)
#plt.subplot(211)
#plt.scatter(caus[:,0],caus[:,1],color='r',marker='.',s=1)
#plt.scatter(src_pt[:100,0],src_pt[:100,1],marker='.',s=1)
#plt.scatter(src_pt[:,0],src_pt[:,1],marker='.',s=1)
plt.plot(caus[:,0],caus[:,1],color='r',label='NN=128')
plt.plot(caus2[:,0],caus2[:,1],color='k',label='NN=32')

plt.scatter(-0.000292478, 0.000334511,color='r',marker='x')
plt.scatter( -0.000292483, 0.000334503,color='k',marker='x')

#plt.scatter(-1.53229e-06, -6.62259e-07,color='k',marker='x')
#plt.scatter( -4.92227e-06, 5.27715e-06,color='r',marker='x')
#plt.plot(crit[:,0][mask],crit[:,1][mask],color='k',label='w sub')
#plt.plot(crit[:,0],crit[:,1],color='k',label='no sub')
'''
plt.axis('equal')
#plt.ylim(-6.0*1e-6,-5.0*1e-6)
#plt.xlim(-4.0*1e-6,-2.0*1e-6)

plt.ylim(-1.0*1e-5,4.0*1e-6)
plt.xlim(-4.0*1e-6,1e-6)
plt.title("NN=64 w sub (seperate)")

plt.subplot(212)
#plt.plot(crit2[:,0][mask2],crit2[:,1][mask2],color='k',label='no sub')

plt.axis('equal')
plt.ylim(-6.0*1e-6,-5.0*1e-6)
plt.xlim(-4.0*1e-6,-2.0*1e-6)
plt.title("NN=64 w/o sub")
'''
plt.legend()
#plt.show()

#
plt.savefig('../../data/glamer/NN_caustic.png')