import os
import numpy as np
import emcee
import triangle as tri

def call_findimage(paras):
    
    ## create input file
        head=open('input.head','r')
        tail=open('input.tail','r')
        inp=open('emcee.input','w')
        
        # read head of input template
        for line in head.readlines():
            #print line
            inp.write(line)
    
        head.close()
        inp.write('\n')
        
        # prepare for lens model parameters
        lens_para=['','','']
        
        for i in range(5):
            lens_para[0]=lens_para[0]+'%f '%paras[i]    #bulge
            lens_para[1]=lens_para[1]+'%f '%paras[5+i]  #disk
                
            if i<2:
                    lens_para[2]=lens_para[2]+'%f '%paras[10+i]  #source
        
        # write lens models into input file
        inp.write('lens sie '+lens_para[0]+'0.0 0.0 \n')
        inp.write('lens sie '+lens_para[1]+'0.0 0.0 \n')
        inp.write('point 1.5 '+lens_para[2]+'\n')
        
        # read & write tail of input template
        for line in tail.readlines():
            inp.write(line)
        
        tail.close()
        inp.close()
        
        ##
        
        ## run glafic
        
        #os.system('glafic')
        #os.system('glafic emcee.input') #this line works for linux machine
        os.system('~/Documents/glafic/glafic emcee.input')
        
        #sub.call('glafic emcee.input')
        
        # read findimg result
        r=np.loadtxt('emcee_point.dat')
        wash=open('emcee_point.dat','w')
        wash.write('0 0 0 0\n 0 0 0 0')
        wash.close()  # wash out after load findimg result

        #print r

        if r[0,0]==4:  # the n_img == 4
            r=r[1:,:] #discard the first line
            r[:,2]=np.absolute(r[:,2]) # ignore parity in magnification
            
            out=np.zeros((len(r[0,:]),3))
                
            for i in range(len(r[:,0])):
                out[i]=np.array(r[i,:3])
        else:
            out=np.array([1]) # flag as not quad system
    


        #print out
        return out   #return findimg result


#########

def img_sort(model):

# sort by distance to x axis/distance to y axis of each image
# model is a 4x3 matrix

        rank=np.zeros(12)
        
        # distance
        
        dis=np.zeros(4)
        for i in range(4):
            dis[i]=model[i,0]**2+model[i,1]**2


        # pick out img D
        d=np.argmin(model[:,1]) # y coord
        rank[9:]=model[d,:]
        model[d,:]=[1000,1000,1000]  # get rid of the recorded set
        

        # sort distance to find image A->B->C
        
        for i in range(3):
            k=np.argmin(dis)
            print model[k,:]
            rank[3*i:3+3*i]=model[k,:]
            dis[k]=1000 #

        return rank


##########


def lnprob(paras):
    
        model=call_findimage(paras)
        
        # B1555 constraints
        #x= np.arrange(10)
        y= np.array([0.0,0.0,17.0,0.0726,0.0480,10.54,0.4117,-0.0280,8.619,0.1619,-0.3680,1.462])
        err= np.array([1.0e-6,1.0e-6,1.70,0.01,0.01,0.97,0.01,0.01,0.83,0.06,0.06,0.13])
        
        # check the n_img
        if len(model) != 1:
            model=img_sort(model)
        
            chi2=0
            for i in range(len(y)):
                chi2=chi2+(model[i]-y[i])**2/err[i]**2
        


        else:
            chi2=1e12

        return -0.5*np.absolute(chi2)


###########
#p=[150.62905,0.17383945,-0.2428983,0.22057935,91.511635,133.6333,0.15735815,-0.2346221,0.8530964,6.221847,0.17383945,-0.2428983]
#p=[150.62905,0.17383945,-0.2428983,0.22057935,91.511635,133.6333,0.15735815,-0.2346221,0.8530964,6.221847,0.0,0.0]

#call_findimage(p)
#print lnprob(p)

nwalker=24 #number of chains
ndim=12   #number of parameters
burn=1000   #number of burn in step
nstep=10000  #number of MCMC steps

## set up for start points

# 2-SIE
mean=[150.62905,0.17383945,-0.2428983,0.22057935,91.511635,133.6333,0.15735815,-0.2346221,0.8530964,6.221847,1.909523e-01,-1.469666e-01]
sig=[1.5,0.005,0.005,0.02,1,2,0.005,0.005,0.005,0.5,0.005,0.005]

w_factor=1.0 # range of random start point

p0=np.zeros((nwalker,ndim))

# random generate star points
for i in range(nwalker):

    c=np.random.rand(ndim)
    
    for j in range(ndim):
        p0[i,j]=mean[j]-w_factor*sig[j]+2*w_factor*sig[j]*c[j]

## sampling

sampler=emcee.EnsembleSampler(nwalker,ndim,lnprob)

## brun-in steps

pos, prob, state = sampler.run_mcmc(p0,burn)

#print pos, prob

## check for valid chain

while 
for i in range(nwalker):
    if prob[i]==-5e11:
        pos[i,j]=mean[j]-w_factor*sig[j]+2*w_factor*sig[j]*c[j] # re-assign starting point

#reset
sampler.reset()

## MCMC run

# how to save each steps--???

#pos, prob, state = sampler.run_mcmc(pos,nstep)


#for i in range(ndim):
#    sampler.flatchain[:,i]
