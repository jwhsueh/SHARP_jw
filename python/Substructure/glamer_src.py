import numpy as np
import matplotlib.pyplot as plt

caus=np.loadtxt('/Users/jwhsueh/Documents/glamer/examples/ParticleExample3/build/caustic.txt')
src_pt=np.loadtxt('/Users/jwhsueh/Documents/glamer/examples/ParticleExample3/build/src_position.txt')

plt.scatter(src_pt[:100,0],src_pt[:100,1],marker='.',s=1)
plt.plot(caus[:,0],caus[:,1],'r')
plt.xlim(-1.0*1e-7,1e-6)
plt.ylim(-1.0*1e-6,5.0*1e-7)
plt.title("N_src=100")
plt.savefig('../../data/glamer/nsrc_100.png')