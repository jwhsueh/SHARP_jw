import os
import numpy as np
import emcee
import itertools as it

img_obs=np.array([0.0,0.0,15.6,-0.075,-0.16,13.104, -1.185,-0.67,6.5208,-1.71,0.460,1.2792])
err_obs=np.array([1.0e-6,1.0e-6,3.70,0.01,0.01,2.97,0.01,0.01,1.83,0.06,0.06,0.3])

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
	    if i<4:
            	lens_para[0]=lens_para[0]+'%f '%paras[i]    #SIE
            	lens_para[1]=lens_para[1]+'%f '%paras[5+i]  #shear

	    else:
		lens_para[0]=lens_para[0]+'%f '%paras[i]    #SIE 5th element
                
            if i<2:
                    lens_para[2]=lens_para[2]+'%f '%paras[9+i]  #source
        
	#print paras
	#print lens_para

        # write lens models into input file
        inp.write('lens sie '+lens_para[0]+'0.0 0.0 \n')
        inp.write('lens pert 1.3390 '+lens_para[1]+'0.0 0.0 \n')
        inp.write('point 1.3390 '+lens_para[2]+'\n')
        
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


        rr=np.array(r[0])
	mock=np.array([0])
	mock=np.append(mock,rr)
        print rr
        print mock
#        print r
	
        nimg=mock[1]

        if nimg==4:  # the n_img == 4
            r=r[1:,:] #discard the first line
            r[:,2]=np.absolute(r[:,2]) # ignore parity in magnification
            
            out=np.zeros((len(r[0,:]),3))
                
            for i in range(len(r[:,0])):
                out[i]=np.array(r[i,:3])
        else:
            out=np.array([1]) # flag as not quad system
    

	wash=open('emcee_point.dat','w')
        wash.write('0 0 0 0\n 0 0 0 0')
        wash.close()  # wash out after load findimg resul

        #print out
        return out   #return findimg result


#########

def img_sort(model):

# sort by distance to x axis/distance to y axis of each image
# model is a 4x3 matrix

        rank=np.zeros(12)
        
        
        # distance
        
        dis=np.zeros((4,4)) # distance to image A,B,C,D
        
        for i in range(4): # model_images
            for j in range(4): # obs_images
                dis[i,j]=(model[i,0]-img_obs[j*3])**2+(model[i,1]-img_obs[j*3+1])**2

        '''
        # sort distance to find image A->B->C->D
        
        for j in range(4):
            k=np.argmin(dis[:,j])
            print model[k,:]
            rank[3*j:3+3*j]=model[k,:]
            dis[k,:]=1000 # kick out from the comparison
        '''

        ## 4x4 combinations, picking up smallest chi^2

        # permutations (24x4 for quad system)
        seq=list(it.permutations([0,1,2,3],4))
        seq=np.array(seq)

        seq_chi=np.zeros(seq.shape[0])

        for i in range(seq.shape[0]):
            for j in range(seq.shape[1]): # image letter (observation)
            # img_obs - img_model in permutation order
                pick=seq[i]
                seq_chi[i]=seq_chi[i]+dis[pick[j],j]

        k=np.argmin(seq_chi)
        
        for i in range(seq.shape[1]):
            pick=seq[k]
            rank[3*i:3*i+3]=model[pick[i],:]

        return rank


##########


def lnprob(paras):
    
        model=call_findimage(paras)

	        

        # B1555 constraints
        #x= np.arrange(10)
        y= img_obs
        err= err_obs
        
        # check the n_img
        if len(model) != 1:
            model=img_sort(model)
		
	    # write findimg chain	    
            writeimg(model)
        
            chi2=0
            for i in range(len(y)):
                chi2=chi2+(model[i]-y[i])**2/err[i]**2
        


        else:
            chi2=np.inf

        return -0.5*np.absolute(chi2)


###########

def start_point(mean,sig,w,nw):
# generate random start point
	
	ndim=len(mean)
    
	pos=np.zeros((nwalker,ndim))
	
	for i in range(nw):
		c=np.random.rand(ndim)
    		for j in range(ndim):
        		pos[i,j]=mean[j]-w_factor*sig[j]+2*w_factor*sig[j]*c[j]

	return pos

###########

def writeimg(paras):
# wrtie the result from findimg into a file

	st=''
	for i in range(len(paras)):
		st=st+'%f '%paras[i]

	find.write(st+'\n')

###########
#p=[150.62905,0.17383945,-0.2428983,0.22057935,91.511635,133.6333,0.15735815,-0.2346221,0.8530964,6.221847,0.17383945,-0.2428983]
#p=[150.62905,0.17383945,-0.2428983,0.22057935,91.511635,133.6333,0.15735815,-0.2346221,0.8530964,6.221847,0.0,0.0]

#call_findimage(p)
#print lnprob(p)

#-------------main-----------------

nwalker=200 #number of chains
ndim=11   #number of parameters
burn=10   #number of (burn in) collect step
burn2=300 #number of (real) burn in step
nstep=5000  #number of MCMC steps

## set up for start points

# SIE+shear
mean=[2.244092e+02,-8.169302e-01,-4.527316e-02,0.4,-6.346850e+01,
-5.842317e-01,  1.003198e-02,  2.491543e-01,  2.265247e+01
,-8.356879e-01,  1.349689e-01]
sig=[15,0.5,0.5,0.2,20,0.5,0.5,0.5,20,0.5,0.5]

w_factor=1.0 # range of random start point


# random generate star points
p0=np.zeros((nwalker,ndim))

p0=start_point(mean,sig,w_factor,nwalker)

## sampling

sampler=emcee.EnsembleSampler(nwalker,ndim,lnprob)


## create findimg chain

find=open('imgchain.dat','w')



## collect valid start point

flag=-np.inf # flag for failure chain

p1=np.zeros(nwalker*ndim)
k=0


'''
#while k<nwalker:
while k<(nwalker-1):
	p0=start_point(mean,sig,w_factor,nwalker) # re-assign starting point
	pos, prob, state = sampler.run_mcmc(p0,burn)
	#print prob
	for i in range(nwalker):
        if prob[i]!=flag:
#			print pos[i,:]			
            p1[(k*ndim):((k+1)*ndim)]=pos[i,:]
            if k==(nwalker-1):
                break
            else:
                k=k+1
    	
	sampler.reset()
	

# reshape p1->p2
p2=np.zeros(((nwalker,ndim)))
for i in range(nwalker):
	p2[i,:]=p1[(i*ndim):((i+1)*ndim)]

# the end of selection step
find.write('# \n')

'''
## (real) burn-in steps

#p3, prob, state = sampler.run_mcmc(p2,burn2)
p3, prob, state = sampler.run_mcmc(p0,burn2)


# the end of burn-in step
find.write('# \n')

## MCMC run
sampler.reset()
#pos, prob, state = sampler.run_mcmc(p2,nstep)

f = open("chain.dat", "w")
g = open('lnprob.dat', 'w')

for result in sampler.sample(p3, iterations=nstep, storechain=False):
    position = result[0]
    lnpro = result[1]
    
    for k in range(nwalker):
	st=''
	st2='%f'%lnpro[k]
	for i in range(ndim):
        	st=st+'%f '%position[k,i]

    	f.write(st+'\n')
	g.write(st2+'\n')
	

f.close()
g.close()
find.close()
