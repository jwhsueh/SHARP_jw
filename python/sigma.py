import numpy as np

file_name=['../data/exp_chain.dat']

k=0
for i in file_name:
    temp=np.loadtxt(i)
    if k==0:
        t=temp
    
    else:
        t=np.append(t,temp,axis=0)
    
    k=k+1


data=t

npara=t.shape[1]
nch=t.shape[0]

mid=nch/2 #median
sigL=np.floor(0.16*nch)
sigR=np.floor(0.84*nch)

for i in range(npara):
    array=t[:,i]
    array=np.sort(array)
    #    print array[sigL],array[mid],array[sigR]
    print array[mid]-array[sigL],array[sigR]-array[mid]

