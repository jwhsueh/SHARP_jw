import numpy as np
import matplotlib.pyplot as plt

path = '/Volumes/sting_1/subs/mock_test/'

smooth_model = np.abs(np.loadtxt(path+'m2_raytrace.txt')[:,8:12])
mock_data = np.abs(np.loadtxt(path+'m1_raytrace.txt')[:,8:12])

drop = np.loadtxt(path+'m2_dropout.txt')
idx = np.arange(0,100)
mask = ~np.in1d(idx,drop)
idx = idx[mask]
print mask[:10]

Rcusp_s = (-smooth_model[:,0]+smooth_model[:,1]-smooth_model[:,2])/(np.abs(smooth_model[:,0])+np.abs(smooth_model[:,1])+np.abs(smooth_model[:,2]))
Rfold_s = (smooth_model[:,1]-smooth_model[:,2])/(np.abs(smooth_model[:,1])+np.abs(smooth_model[:,2]))

Rcusp_s,Rfold_s = Rcusp_s[mask],Rfold_s[mask]
print Rcusp_s[:10]

Rcusp_m = (-mock_data[:,0]+mock_data[:,1]-mock_data[:,2])/(np.abs(mock_data[:,0])+np.abs(mock_data[:,1])+np.abs(mock_data[:,2]))
Rfold_m = (-mock_data[:,1]+mock_data[:,2])/(np.abs(mock_data[:,1])+np.abs(mock_data[:,2]))
Rcusp_m,Rfold_m = Rcusp_m[mask],Rfold_m[mask]

plt.hist(Rcusp_s)
plt.hist(Rcusp_m,color='r',alpha=0.6)
#plt.scatter(idx,np.abs(Rcusp_s))
#plt.scatter(idx,np.abs(Rcusp_m),color='r')
plt.show()