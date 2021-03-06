import os
import numpy as np
import emcee
import scipy.stats
import itertools as it
import commands

img_obs=np.array([0.0,0.0,17.0,-0.0726,0.0480,10.54,-0.4117,-0.0280,8.619,-0.1619,-0.3680,1.462])
err_obs=np.array([1.0e-6,1.0e-6,1.70,0.01,0.01,0.97,0.01,0.01,0.83,0.06,0.06,0.13])


def call_findimage(paras):
    
    ## create input file
        head=open('input_exp.head','r')
        #tail=open('input_exp.tail','r')
        inp=open('B1555_emcee.input','w')
        
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
            #lens_para[1]=lens_para[1]+'%f '%paras[5+i]  #disk


        lens_para[1]=lens_para[1]+'%f '%paras[5]+'-1.615516e-1 -2.056755e-1 8.28e-1 -8.0e+0 0.0 0.0 '+'%f '%paras[6] #r_s for expdisk
        
        # write lens models into input file
        inp.write('alpha '+lens_para[0]+'0.0 0.0 0.0 0.0 1.0 \n') #SIE
        inp.write('expdisk '+lens_para[1]+'0.0 0.0 \n')       #expdisk

        inp.write('\n 0 0 0 0 0 0 0 0 0 0 \n 0 0 0 0 0 0 0 0 0 0 \n')

        # findimg [souce pos]

        lens_para[2]=lens_para[2]+'%f '%paras[7]+'%f '%paras[8]

        inp.write('\n findimg '+lens_para[2]+'\n')
        

        inp.close()
    
        print_out=commands.getstatusoutput('./lensmodel B1555_emcee.input')
        print_out=print_out[1]
        print_out=print_out.split('\n')
                        
        ot=len(print_out)

        # extract findimg result

        # get img number info

        check=print_out[-6]
        #print check
        check=check.split()
        #print check[0]

        if check[0] == '#':

            f_str=[]
            out=np.zeros(12)
            j=0
            for i in [5,4,3,2]:

                f_str.append(print_out[-i])
                sp=f_str[-1].split()
                out[3*j:3*j+3]=[np.float(sp[0]),np.float(sp[1]),np.float(sp[2])]
                j=j+1

#           for i in range(4):
#                out[3*i:3*i+3]=[np.float(f_str[4*i]),np.float(f_str[4*i+1]),np.float(f_str[4*i+2])]


        else:
            out=np.array([1])


        #print out
        return out   #return findimg result


#########

def img_sort(model):

# pick the smallest chi^2 as image sorting standard
# model is a 1x12 matrix

    model=model.reshape((4,3))
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

    #print rank
    return rank


##########


def lnprob(paras):
    
    #print paras
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

    print chi2
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

nwalker=200 #number of chains
ndim=9   #number of parameters
burn=10   #number of (burn in) collect step
burn2=1000 #number of (real) burn in step
nstep=10000  #number of MCMC steps

## set up for start points

# expdisk
mean=[1.776548e-01, -1.772157e-01, -2.049533e-01, 2.500000e-01, -4.890300e+00,
      2.734703e-01, 2.355863e-01,
      -1.971670e-01, -1.512136e-01]
#sig=[0.01,0.01,0.01,0.01,1,0.01,0.01,0.01,0.01,1,0.01,0.01,0.01]
sig=[0.02,0.01,0.01,0.2,5,0.1,0.02,0.01,0.01]
#sig=np.zeros(13)

w_factor=1.0 # range of random start point


# random generate star points
p0=np.zeros((nwalker,ndim))

p0=start_point(mean,sig,w_factor,nwalker)

## sampling

sampler=emcee.EnsembleSampler(nwalker,ndim,lnprob)


## create findimg chain

find=open('exp_imgchain.dat','w')


## brun-in steps

#pos, prob, state = sampler.run_mcmc(p0,burn)


## (real) burn-in steps

#p3, prob, state = sampler.run_mcmc(p2,burn2)
p3, prob, state = sampler.run_mcmc(p0,burn2)


# the end of burn-in step
find.write('# \n')


## MCMC run
sampler.reset()
#pos, prob, state = sampler.run_mcmc(p2,nstep)

f = open("exp_chain.dat", "w")
g = open('exp_lnprob.dat', 'w')

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

