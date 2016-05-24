import numpy as np
import matplotlib.pyplot as plt

AngMomIndex = np.loadtxt('../../data/illustris_1/kinematics_99.dat',dtype = 'float',unpack=True, usecols=[3])

plt.hist(AngMomIndex,50)
plt.show()