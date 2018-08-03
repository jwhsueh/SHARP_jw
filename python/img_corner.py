import numpy as np
import corner
import matplotlib.pyplot as plt
#import matplotlib.patches as mpat

fsub='0500'
lens = 'mock3'

#path = '../data/sub_gravlens/'
#t = np.loadtxt(path+'B1422_flatchain_0.txt')
#t = np.loadtxt(path+'B1422_eta.txt')
#w = np.loadtxt(path+'mcmc_chi2_1.txt')[:50000]

#t = np.loadtxt('/Volumes/sting_1/subs/MG0414/result_new/MG0414_0100_1out.txt')
t = np.loadtxt('/Volumes/sting_1/subs/'+lens+'/result_new3/'+lens+'_'+fsub+'_out_com.txt')
#chi2 = np.loadtxt('/Volumes/sting_1/subs/result/B1422_sub01_per_chi2.txt')
#w = np.exp(-1.0*chi2/2.0)
#w = w/np.sum(w)
#t = np.loadtxt('/Volumes/sting_1/subs/glamer_out/copy.txt')

f_idx = 8+2


print t.shape

'''
table = open(path+'B1422_realization1/smooth_result.txt','r')
lines = table.readlines()
#print lines[0].split()
i = 0
t = np.empty((5000,9))
for Aline in lines:
    arr = np.array(Aline.split()).astype(float)
    t[i,:] = arr
    i = i+1

print t[:10,:]
'''
'''
k=0
for i in file_name:
    temp=np.loadtxt(i)
    if k==0:
        t=temp

    else:
        t=np.append(t,temp,axis=0)

    k=k+1
'''
'''
xa,ya,fa=t[:,0],t[:,1],t[:,2]
xb,yb,fb=t[:,3],t[:,4],t[:,5]
xc,yc,fc=t[:,6],t[:,7],t[:,8]
xd,yd,fd=t[:,9],t[:,10],t[:,11]

rab=fa/fb
rac=fa/fc
rad=fa/fd
'''
#data=t[:,:-1]
data=t
data[:,8:12]=np.abs(data[:,8:12])
#f = data[:,f_idx]
data[:,8]=data[:,8]/data[:,f_idx]
data[:,9]=data[:,9]/data[:,f_idx]
#data[:,10]=data[:,10]/data[:,f_idx]
data[:,11]=data[:,11]/data[:,f_idx]

data = np.delete(data,f_idx,1)
data=data[:,:11]
print data[0,:]
## triangle plot

#ndim, nsamples = 10, len(chi2)

figure = corner.corner(data,
   #labels=[r"$b$", r"$xc$",r"$yc$", r"$e$",r"$PA$",r"gamma1",r"$gamma2$", r"$sx$",r"$sy$"],
   #labels=[r"$xa$", r"$ya$",r"$fb/fa$",r"xc",r"$yc$", r"$fc/fa$",r"$xd$", r"$yd$",r"$fd/fa$"],#weights=w,
   range=[(0.9499,0.951),(-0.01,0.01),(0.38,0.39),(-0.35,-0.33),(-0.805,-0.80),(-0.01,0.01),(0.31,0.325),(-0.75,-0.74),(0.03,0.05),(1.0,1.23),(0.5,0.7)],
   truths = [0.95065,0.0,0.38925,-0.33388,-0.80215,0.0,0.31998,-0.74771,0.03412,1.04549,0.556], #mock3
   #quantiles=[0.16, 0.5, 0.84],
   show_titles=True, verbose=True)
   #truths = [7.466419e-01, 7.645039e-01, -6.730659e-01, 3.782159e-01 ,5.556257e+01 ,1.473030e-01, 5.249074e+01,3.902599e-01, -4.156310e-01])
   #truths = [7.466419e-01,7.645039e-01, -6.730659e-01, 3.782159e-01 ,5.556257e+01, 1.473030e-01 ,5.249074e+01,3.902599e-01, -4.156310e-01])
   #truths = [7.411447e-01,7.655555e-01,-6.742337e-01,3.793138e-01,5.560879e+01,1.473229e-01,5.272611e+01,3.870229e-01,-4.128895e-01])
   #truths = [0.7419989,7.703456e-1,-6.776905e-1,3.781747e-1,55.64814,1.487065e-1,52.77599,3.851740e-1,-4.127741e-1])
   #weights=none)
    #range=[(0.75,0.8),(-0.7,-0.68),(-0.7,-0.65),(0.35,0.45),(52.,57),(0.11,0.16),(51.,54.),(0.37,0.4),(-0.41,-0.4)])

#                         quantiles=[0.16, 0.5, 0.84],
#                         show_titles=True, title_args={"fontsize": 12})
                         
#figure.savefig("../data/sub_gravlens/B1422_sam_corner.png")
figure.savefig("/Volumes/sting_1/subs/"+lens+"/result_new3/"+lens+"_"+fsub+"_com_out.png")





