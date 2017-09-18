import numpy as np
import gravlens_tool as gt
import RcuspTool as rcusp
import matplotlib.pyplot as plt

subname = np.genfromtxt('/Volumes/sting_1/data/snap99_elp2_Rcusp.txt',dtype='str')

table = np.empty((len(subname),4))
for i in range(len(subname)):
	filename = '../../data/sub_gravlens/snap99_elp/'+subname[i]+'_best.dat'

	img_x,img_y,img_f = gt.get_bestdat(filename)
	#print img_x,img_y,img_f

	table[i,:] = rcusp.rcusp_tool(img_x,img_y,img_f)
	#print table

#print table
plt.scatter(table[:,2],np.abs(table[:,1]))
plt.show()