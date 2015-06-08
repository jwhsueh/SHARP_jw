import os
import numpy as np
import emcee
import itertools as it
import commands

img_obs=np.array([0.0,0.0,15.6,-0.075,-0.16,13.104, -1.185,-0.67,6.5208,-1.71,0.460,1.2792])
err_obs=np.array([1.0e-6,1.0e-6,3.70,0.01,0.01,2.97,0.01,0.01,1.83,0.06,0.06,0.3])

# center of the images
xc=(img_obs[0]+img_obs[3]+img_obs[6]+img_obs[9])/4.0
yc=(img_obs[1]+img_obs[4]+img_obs[7]+img_obs[10])/4.0

b=0.5 # limitation radius of component position



def call_calcimage(paras):
    
    #print paras
    
    ## check on paras physical limitation
    # send unqualify model to high chi^2
    
    # ellipticity
    if paras[3]>1.0 or paras[3]<0.0:
        out=np.array([1]) # flag as unqualify system
    
    elif paras[7]>1.0 or paras[7]<0.0:
        out=np.array([1])

    # velocity dispersion
    elif paras[0]<0.0:
        out=np.array([1])
    
    # 1st component position
    elif (paras[1]-xc)**2+(paras[2]-yc)**2>b**2:
        out=np.array([1])

    # 2nd component position
    elif (paras[5]-xc)**2+(paras[6]-yc)**2>b**2:
        out=np.array([1])
        
    # source position
    elif (paras[9]-xc)**2+(paras[10]-yc)**2>b**2:
        out=np.array([1])
    
    else:
        out=run_glafic(paras)



    return out   #return calcimage result (source positions from each obs img)


##########

def run_glafic(paras):
    ## create input file
    head=open('input_src.head','r')
    tail=open('input_src.tail','r')
    inp=open('emcee_src.input','w')
        
    # read head of input template
    for line in head.readlines():
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



    # write lens models into input file
    inp.write('lens sie '+lens_para[0]+'0.0 0.0 \n')
    inp.write('lens pert 1.3390 '+lens_para[1]+'0.0 0.0 \n')
    inp.write('point 1.3390 '+lens_para[2]+'\n')
        
    # read & write tail of input template
    for line in tail.readlines():
        inp.write(line)

    tail.close()
    inp.close()
        

        
    ## run glafic
        
    #os.system('glafic')
    #os.system('glafic emcee.input') #this line works for linux machine
    #os.system('~/Documents/glafic/glafic emcee.input')
        
    print_out=commands.getstatusoutput('~/Documents/glafic/glafic emcee_src.input')
    print_out=print_out[1]
    print_out=print_out.split('\n')
        
    ot=len(print_out)
        
    # extract src position [A,B,C,D]
    src_st=[]
    for i in [3,2,1,0]:
        src_st.append(print_out[ot-13*i-2]) # src_x
        src_st.append(print_out[ot-13*i-1]) # src_y
        
    # split string
    out=np.zeros(8)
    for i in range(8):
        
        sp=src_st[i].split()
        print sp
        out[i]=np.float(sp[2])


    return out

#########

'''
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


        # sort distance to find image A->B->C->D
        
        #for j in range(4):
        #    k=np.argmin(dis[:,j])
        #    print model[k,:]
        #    rank[3*j:3+3*j]=model[k,:]
        #    dis[k,:]=1000 # kick out from the comparison
        

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
'''

def lnprob(paras):
    
        model=call_calcimage(paras)
        # model is a 8x1 matrix for src positions
	        

        # src plane constraints
        #x= np.arrange(10)
        #y= img_obs
        #err= err_obs
        
	    # write findimg chain	    
        writeimg(model)
        
        
        # source plane penalty function
        chi2=0
        x=np.zeros(4)
        y=np.zeros(4)
        
        for i in range(4):
            x[i]=model[2*i]
            y[i]=model[2*i+1]
        
        x=np.array(x)
        y=np.array(y)
        x0=np.mean(x)
        y0=np.mean(y)
        
        for i in range(4):
            chi2=chi2+(x[i]-x0)**2+(y[i]-y0)**2
        


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
mean=[2.324253e+02,-1.100387e+00,-9.439666e-02,0.8,60,
-1.068075e-01, -2.923237e-01, 0.1,  9.857114e+00
,-6.550547e-01,  2.884248e-01]
sig=[10,0.2,0.2,0.1,10,0.2,0.2,0.05,10,0.2,0.2]

w_factor=1.0 # range of random start point


# random generate star points
p0=np.zeros((nwalker,ndim))

p0=start_point(mean,sig,w_factor,nwalker)


## sampling

sampler=emcee.EnsembleSampler(nwalker,ndim,lnprob)


## create findimg chain

find=open('srcchain.dat','w')



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

f = open("chain_src.dat", "w")
g = open('lnprob_src.dat', 'w')

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
