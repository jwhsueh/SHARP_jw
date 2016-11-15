import numpy as np
import matplotlib.pyplot as plt

caus=np.loadtxt('/Users/jwhsueh/Documents/glamer/examples/ParticleExample/build/caustic_wsub.txt')
caus2=np.loadtxt('/Users/jwhsueh/Documents/glamer/examples/ParticleExample/build/caustic_nosub.txt')
src_pt=np.loadtxt('/Users/jwhsueh/Documents/glamer/examples/ParticleExample/build/src_position2.txt')

caus_intv=np.sqrt((caus[1:,0]-caus[:-1,0])**2+(caus[1:,1]-caus[:-1,1])**2)
print caus[1:,0]-caus[:-1,0]
print min(caus_intv)

#plt.scatter(caus[:,0],caus[:,1],color='r',marker='.',s=1)
#plt.scatter(src_pt[:100,0],src_pt[:100,1],marker='.',s=1)
#plt.scatter(src_pt[:,0],src_pt[:,1],marker='.',s=1)
plt.plot(caus[:,0],caus[:,1],color='r',label='w/ sub')
plt.plot(caus2[:,0],caus2[:,1],color='k',label='no sub')

#plt.xlim(-1.0*1e-6,1e-6)
#plt.ylim(-1.0*1e-6,5.0*1e-6)
#plt.title("N_src=100")
#plt.show()
plt.legend()
plt.savefig('../../data/glamer/caustics.png')