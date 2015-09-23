import os
import numpy as np
import emcee
import scipy.stats
import itertools as it
import commands

img_obs=np.array([0.0,0.0,17.0,-0.0726,0.0480,10.54,-0.4117,-0.0280,8.619,-0.1619,-0.3680,1.462])
err_obs=np.array([1.0e-6,1.0e-6,1.071,0.001,0.001,0.664,0.001,0.001,0.543,0.006,0.006,0.092])

Er=0.26 #Einstein radius in arcsec
center=np.array([-0.1875,-0.1625]) #Einstein ring center
f_lim=np.array([0.001,0.1]) #upper limit of f_sub


def call_findimage(paras):
    
	paras[-3]=10**(paras[-3]) # b_sub

#print paras
	
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
        lens_para=['','','','']
        
        for i in range(5):
            lens_para[0]=lens_para[0]+'%f '%paras[i]    #bulge
            lens_para[1]=lens_para[1]+'%f '%paras[5+i]  #disk

	for i in [-3,-2,-1]:
		lens_para[2]=lens_para[2]+'%f '%paras[i] #substructure

        
        lens_para[1]=lens_para[1]+'0.0 0.0 ''%f '%paras[10] #r_s for expdisk

	# write real chain

	mcs=''
	for i in range(len(paras)):
	        print paras[i]
		mcs=mcs+'%f '%paras[i]

        mcmc.write(mcs+'\n') # order: comp1+comp2+src+sub
 
        
        # write lens models into input file
        inp.write('alpha '+lens_para[0]+'0.0 0.0 0.0 0.0 1.0 \n') #SIE
        inp.write('expdisk '+lens_para[1]+'0.0 0.0 \n')       #expdisk
	inp.write('alpha '+lens_para[2]+'0.0 0.0 0.0 0.0 0.0 0.0 1.0 \n')

        inp.write('\n 0 0 0 0 0 0 0 0 0 0 \n 0 0 0 0 0 0 0 0 0 0 \n 0 0 0 0 0 0 0 0 0 0 \n')

        # findimg [souce pos]

        lens_para[3]=lens_para[3]+'%f '%paras[11]+'%f '%paras[12]

        inp.write('\n findimg '+lens_para[3]+'\n')
        

        inp.close()
        
        ##
        
        ## run gravlens
        
        #os.system('glafic')
        #os.system('glafic emcee_exp.input') #this line works for linux machine
        #os.system('./lensmodel emcee_exp.input')

        
        # read findimg result [output on screen, gravlens]

        print_out=commands.getstatusoutput('~/Documents/gravlens/lensmodel B1555_emcee.input')
        print_out=print_out[1]
        print_out=print_out.split('\n')
                        
        ot=len(print_out)

        ## extract findimg result

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
    	else:
        	out=np.array([1])

	return out   #return findimg result


#########

def check_paras(paras):
	c=0
	## source pos
	xs,ys=paras[11],paras[12]
	x1,y1=paras[1],paras[2]
	x2,y2=paras[5],paras[6]
	xsb,ysb=paras[14],paras[15]
	bsub=10**paras[13]

	# a list of check criteria
	#cc=np.array([xs**2+ys**2,xsb])
	

	if (xs-center[0])**2+(ys-center[1])**2>Er**2:
		c=1 # flag of invalid source position
		
		
	elif bsub>f_lim[1]*Er: 
		c=1
		
		
	elif (x1-center[0])**2+(y1-center[1])**2>Er**2:
		c=1 # flag of invalid source position
		
	elif (x2-center[0])**2+(y2-center[1])**2>Er**2:
		c=1 # flag of invalid source position

	return c

	


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
    
	## criteria check for paras
	p_flag=check_paras(paras)

	if p_flag !=1:
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

	else:
		chi2=np.inf

    #print chi2
	prob_chi2.write('%f'%chi2+'\n')
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
ndim=16   #number of parameters
burn=10   #number of (burn in) collect step
burn2=1000 #number of (real) burn in step
nstep=10000  #number of MCMC steps

## set up for start points

# expdisk
mean=[2.052429e-01, -1.818271e-01, -1.987580e-01, 3.091056e-01, 2.957686e+00,
      1.100000e-01, -1.471605e-01, -2.056755e-01, 8.660612e-01, 6.996249e+00,2.770794e-01,
      -1.952576e-01, -1.493771e-01,-2.585,center[0],center[1]]
#sig=[0.01,0.01,0.01,0.01,1,0.01,0.01,0.01,0.01,1,0.01,0.01,0.01]
sig=[0.05,0.05,0.05,0.02,10,0.05,0.05,0.05,0.02,10,0.05,0.05,0.05,1.0,Er/2,Er/2]
#sig=np.zeros(13)

w_factor=1.0 # range of random start point


# random generate star points
p0=np.zeros((nwalker,ndim))

p0=start_point(mean,sig,w_factor,nwalker)
#print p0

## sampling

sampler=emcee.EnsembleSampler(nwalker,ndim,lnprob)


## create findimg chain

find=open('B1555_imgchain.dat','w')

## create real para chain & prob chain
mcmc=open('B1555_parachain.dat','w')
prob_chi2=open('B1555_chi2chain.dat','w')

## brun-in steps

#p3, prob, state = sampler.run_mcmc(p2,burn2)
p3, prob, state = sampler.run_mcmc(p0,burn2)


# the end of burn-in step
find.write('# \n')
mcmc.write('# \n')
prob_chi2.write('# \n')

# I'm try to flush it
find.close()
mcmc.close()
prob_chi2.close()

## create new write-in file again

find=open('B1555_imgchain.dat','w')

mcmc=open('B1555_parachain.dat','w')
prob_chi2=open('B1555_chi2chain.dat','w')

## MCMC run
sampler.reset()
#pos, prob, state = sampler.run_mcmc(p2,nstep)

for result in sampler.sample(p3, iterations=nstep, storechain=False):
    position = result[0]
    lnpro = result[1]
	

find.close()
mcmc.close()
prob_chi2.close()
