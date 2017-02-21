import numpy as np
import plot_lensmod as pltlm
import matplotlib.pyplot as plt

lens_path = '/Users/jwhsueh/Documents/SHARP_jw/models/snap99_179899/'

crit_file = lens_path+'crit.dat'

## read-in magnification

mag_3 = np.zeros(280)
sx, sy = np.zeros(280),np.zeros(280)

for i in range(280):
	table = np.loadtxt(lens_path+'/rt_gravlens/rt_src'+str(i)+'.txt',skiprows=1)

	mag = table[:,2]
	mag = list(mag)
	mag.pop(mag.index(np.min(mag)))
	mag_3[i] = np.min(mag)

	table = open(lens_path+'/rt_gravlens/rt_src'+str(i)+'.txt','r')
	line = table.readline().split()
	sx[i] = float(line[0])
	sy[i] = float(line[1])
	#print sx[i],sy[i]


plt.contour(sx,sy,mag_3)
plt.show()
