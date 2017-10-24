import numpy as np
import matplotlib.pyplot as plt

#real_list = np.array([138,  589,  915, 1460, 3703, 3841, 3864, 4882,5106, 5461, 5768, 5928, 6077, 6862, 8023, 8850]).astype(str)
real_list = np.array([0]).astype(str)
for  realid in real_list:
	los_table = np.loadtxt('/Volumes/sting_1/subs/real_01/real'+realid+'.txt')
	#z=los_table[:,0]
	mass = los_table[:,0]
	x = los_table[:,1]+7.65020767e-01
	y = los_table[:,2]-6.73284706e-01

	img_x = np.array([0.0,0.38925,-0.33388,0.95065])
	img_y = np.array([0.0,0.31998,-0.74771,-0.80215])

	plt.scatter(img_x,img_y,marker='*',color='r',s=50)
	plt.scatter(x,y,marker='o',facecolor='none',edgecolor='b',s=mass/1e6)
	plt.title('chi2<10')
	plt.gca().set_aspect('equal')
	plt.show()
	#plt.savefig('/Volumes/sting_1/subs/real_01/real'+realid+'.png')
	#plt.clf()