import numpy as np
import matplotlib.pyplot as plt
import gravlens_tool as gt
import emcee
import sys

#case_idx = int(sys.argv[1])
lens = 'mock2'

path = '../../data/sub_gravlens/'
#path = '../data/sub_gravlens/'
imgfile = lens+'_findimg.input' #findimg

#table = np.loadtxt('../../data/sub_gravlens/B1422_realization1/real0.txt')
#table = np.loadtxt('../data/sub_gravlens/B1422_realization1/real0.txt')
#x_list,y_list,bsub_list,rt_list = table[:,1],table[:,2],table[:,3],table[:,4]

zl_g = 0.6 
zs_g = 3.12
zg = np.array([zl_g,zs_g])

## model parameter
#macro_mod = np.array([[1.118581e+00, -5.571237e-01, -1.349607e+00, 3.976823e-01, -7.274792e+01,3.629423e-02, 7.004755e+01],
    #[5.880343e-02, -1.980000e-01, 1.800000e-01, 0.0, 0.0, 0.0, 0.0]])

macro_mod = np.array([[7.537347e-01, 7.457579e-01, -6.602210e-01, 3.239710e-01, 5.596422e+01 ,1.616563e-01, 5.258606e+01],
    [1e-11, 2.0,2.0, 0.0, 0.0, 0.0, 0.0]])

#print macro_mod[0,:]

src = np.array([3.887401e-01, -4.143577e-01])


obs_x = np.array([9.507944e-01,  1.493761e-03,   3.892554e-01, -3.340559e-01])
obs_y = np.array([-8.022173e-01 , 1.813538e-03,  3.199820e-01 ,-7.484454e-01])
obs_f = np.array([0.03412,1.198,1.04549,0.556])
idx_f = 2


err_x = np.array([0.003,0.003,0.003,0.003])
err_y = np.array([0.003,0.003,0.003,0.003])
err_f = obs_f*0.05

burn_flag = 0

findimg_outfile = open(path+lens+'_mcmc_mod.txt','w')
chi2_outfile = open(path+lens+'_mcmc_chi2.txt','w')
param_outfile = open(path+lens+'_mcmc_param_lens.txt','w')

### --- findimg process

def do_findimg(paras):
    if np.abs(paras[0])<2.0:

        step_mod_g = np.array([np.abs(paras[0]),paras[1],paras[2],paras[3],paras[4],paras[5],paras[6]])
        step_mod = np.array([step_mod_g,macro_mod[1,:]]) ## edit here when there're luminous satellite
        #step_mod = step_mod_g

        step_src = np.array([paras[7],paras[8]])

        gt.create_findimg_macro(step_mod,['sie','sie'],zg,step_src,path,imgfile)
        qflag = gt.run_findimg(path,lens,imgfile) # flag tells you if it's quad

        if qflag==True:
            mod_x,mod_y,mod_f = gt.get_imgresult(path,lens)
            #print mod_x
            #mod_x,mod_y,mod_f = findimg_sort(mod_x,mod_y,mod_f)
            mod_f = np.abs(mod_f)
            mod_f = mod_f/mod_f[idx_f]
            #print mod_f,obs_f
                
            chi2=0
            for i in range(4):
                chi2=chi2+(mod_x[i]-obs_x[i])**2/err_x[i]**2+(mod_y[i]-obs_y[i])**2/err_y[i]**2+(mod_f[i]-obs_f[i])**2/err_f[i]**2

            chi2 = chi2/9.0

            if (burn_flag!=0):
                
                if chi2<=60.:
                    print mod_x,mod_y,mod_f

                    chi2_outfile.write(str(chi2)+'\n')
                    param_outfile.write(str(np.abs(paras[0]))+'\t'+str(paras[1])+'\t'+str(paras[2])+'\t'+str(paras[3])+'\t'+str(paras[4])+'\t'+str(paras[5])+'\t'+str(paras[6])
                        +'\t'+str(paras[7])+'\t'+str(paras[8])+'\n')
        else:
            chi2=100000

    else:
        chi2=100000

    print 'chi2 = '+str(chi2)

    return chi2


### --- sort the findimg output

def findimg_sort(mod_x,mod_y,mod_f):
    xsort = np.sort(mod_x)

    x_new,y_new,f_new = np.empty(4),np.empty(4),np.empty(4)
    for i in range(4):
        idx = list(mod_x).index(xsort[i])

        x_new[i],y_new[i],f_new[i] = mod_x[idx],mod_y[idx],mod_f[idx]

    return x_new,y_new,f_new

### --- log prob for emcee
def lnprob(paras):
    chi2 = do_findimg(paras)
    return -0.5*np.absolute(chi2)


### --- uniform random starting point

def start_point(mean,sig,w,nw):
# generate random start point
	
	ndim=len(mean)
    
	pos=np.zeros((nwalker,ndim))
	
	for i in range(nw):
		c=np.random.rand(ndim)
    		for j in range(ndim):
        		pos[i,j]=mean[j]-w_factor*sig[j]+2*w_factor*sig[j]*c[j]

	return pos


## --- emcee run settings

nwalker=50 #number of chains
ndim=9   #number of parameters
burn=200   #number of (burn in) collect step
nstep=1000  #number of MCMC steps

## set up for start points

#mean=[macro_mod[0,0], macro_mod[0,1],macro_mod[0,2],macro_mod[0,3],macro_mod[0,4],macro_mod[0,5],macro_mod[0,6],src[0],src[1]]
mean=[macro_mod[0,0], macro_mod[0,1],macro_mod[0,2],macro_mod[0,3],macro_mod[0,4],macro_mod[0,5],macro_mod[0,6],src[0],src[1]]
#sig=[0.01,0.01,0.01,0.01,1,0.01,0.01,0.01,0.01,1,0.01,0.01,0.01]
sig=np.array([0.001,0.001,0.01,0.01,1,0.001,1,0.001,0.001])
#sig=np.zeros(13)

w_factor=1.0 # range of random start point


# random generate star points
#p0=np.zeros((nwalker,ndim))

p0=start_point(mean,sig,w_factor,nwalker)

## sampling

sampler=emcee.EnsembleSampler(nwalker,ndim,lnprob)

## burn-in stage
print '### burn-in stage started ###'
print '\n'

pos, prob, state = sampler.run_mcmc(p0,burn)

burn_flag = 1
print '### burn-in stage finished ###'
print '\n'
print '### start MCMC run ###'

## MCMC run
sampler.reset()
sampler.run_mcmc(pos,nstep)

chi2_outfile.close()
param_outfile.close()
'''
f = open("exp_chain.dat", "w")
g = open('exp_lnprob.dat', 'w')

for result in sampler.sample(p0, iterations=nstep, storechain=False):
    print 'in the loop'
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
'''


## --- save the flat chain
np.savetxt(path+lens+'_flatchain.txt',sampler.flatchain)

