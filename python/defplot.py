import matplotlib.pyplot as plt
import numpy as np

t=np.loadtxt('./gravlens/def_2.dat')

#print t[:,0]
fig=plt.quiver(-1.0*t[:,0],t[:,1],t[:,3],t[:,4])
plt.show()
