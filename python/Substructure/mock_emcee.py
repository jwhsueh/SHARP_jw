import numpy as np
import matplotlib.pyplot as plt
import gravlens_tool as gt
import emcee
import sys

case_idx = int(sys.argv[1])

path = '../../data/sub_gravlens/'
#path = '../data/sub_gravlens/'
imgfile = 'B1422_real_findimg_'+str(case_idx)+'.input' #findimg

table = np.loadtxt('../../data/sub_gravlens/B1422_realization1/real0.txt')
#table = np.loadtxt('../data/sub_gravlens/B1422_realization1/real0.txt')
x_list,y_list,bsub_list,rt_list = table[:,1],table[:,2],table[:,3],table[:,4]

zl_g = 0.34
zs_g = 3.62

class macro_mod: ## B1422
    b =  7.466390e-01 
    xc =7.645506e-01
    yc =  -6.730964e-01
    e =  3.782699e-01
    PA = 5.556224e+01 
    gamma1 =1.473014e-01
    gamma2 = 5.249086e+01 ## shear PA
    sx =3.902577e-01 # src position
    sy = -4.156271e-01
'''
class macro_mod: ## B1422
    b =  0.7419989
    xc =7.703456e-1
    yc =  -6.776905e-1
    e =  3.781747e-1
    PA = 55.64814
    gamma1 =1.487065e-1
    gamma2 = 52.77599## shear PA
    sx =3.851740e-1 # src position
    sy = -4.127741e-1
'''
'''
class micro_mod:
	bi = bsub_list
	rt_i = rt_list
	xi = x_list
	yi = y_list
'''
class micro_mod:
    ks_i = np.zeros(1)
    rs_i = np.zeros(1)
    xi = np.zeros(1)
    yi = np.zeros(1)


### ---- here's the part that changes realization
if case_idx==0: ## real observation
    obs_x = np.array([-3.338751e-01,0,3.892583e-01,9.506678e-01])
    obs_y = np.array([-7.477099e-01,0, 3.199784e-01,-8.021654e-01])
    obs_f = np.array([0.551,1.062, 1.,0.024])


elif case_idx==1:
    obs_x = np.array([-0.3339,    0.,0.3893,      0.9507])
    obs_y = np.array([-0.7477,    -0.,0.32,      -0.8022])
    obs_f = np.array([4.313715  ,7.464847,6.038479,2.485514e-01])

elif case_idx==2:
    obs_x = np.array([-0.3339,   0. ,0.3893 ,     0.9507])
    obs_y = np.array([  -0.7477,   -0. ,0.32  ,    -0.80221])
    obs_f = np.array([3.893091,7.352067 ,5.637157,2.561626e-01 ])
    obs_f = obs_f/np.max(obs_f)

elif case_idx==3:
    obs_x = np.array([-0.3339,   0. ,0.3893 ,     0.9507])
    obs_y = np.array([  -0.7477,   -0. ,0.32  ,    -0.80221])
    obs_f = np.array([3.985744,7.416987 ,6.048564,2.467289e-01])
    obs_f = obs_f/np.max(obs_f)

elif case_idx==4:
    obs_x = np.array([-0.3339,   0. ,0.3893 ,     0.9507])
    obs_y = np.array([  -0.7477,   -0. ,0.32  ,    -0.80221])
    obs_f = np.array([3.600924,6.102982,5.030562,2.886554e-01 ])
    obs_f = obs_f/np.max(obs_f)


err_x = np.array([0.00005,0.00005,0.00005,0.00005])
err_y = np.array([0.00005,0.00005,0.00005,0.00005])
err_f = np.array([0.007,0.009,1e-5,0.006])

burn_flag = 0

findimg_outfile = open(path+'mcmc_mod'+str(case_idx)+'.txt','w')
chi2_outfile = open(path+'mcmc_chi2_'+str(case_idx)+'.txt','w')
param_outfile = open(path+'mcmc_param_'+str(case_idx)+'.txt','w')

### --- findimg process

def do_findimg(paras):
    #print 'wtf'
    class step_mod:
        b = paras[0]
        xc = paras[1]
        yc = paras[2]
        e = paras[3]
        PA = paras[4]
        gamma1 = paras[5]
        gamma2 = paras[6]
        sx = paras[7]
        sy = paras[8]
        zl = zl_g
        zs = zs_g
    gt.create_findimg(step_mod,micro_mod,path,imgfile)
    qflag = gt.run_findimg(path,imgfile,case_idx) # flag tells you if it's quad

    if qflag==True:
        mod_x,mod_y,mod_f = gt.get_imgresult(path,case_idx)
        mod_x,mod_y,mod_f = findimg_sort(mod_x,mod_y,mod_f)
        mod_f = np.abs(mod_f)
        mod_f = mod_f/mod_f[2]
        #print mod_f,obs_f
            
        chi2=0
        for i in range(4):
            chi2=chi2+(mod_x[i]-obs_x[i])**2/err_x[i]**2+(mod_y[i]-obs_y[i])**2/err_y[i]**2+(mod_f[i]-obs_f[i])**2/err_f[i]**2

        chi2 = chi2/9.0

        if (burn_flag!=0):
            
            if chi2<=60.:
                print mod_x,mod_y,mod_f

                chi2_outfile.write(str(chi2)+'\n')
                param_outfile.write(str(paras)[1:-1]+'\n')
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

# expdisk
mean=[macro_mod.b, macro_mod.xc,macro_mod.yc,macro_mod.e,macro_mod.PA,macro_mod.gamma1,macro_mod.gamma2,macro_mod.sx,macro_mod.sy]
#sig=[0.01,0.01,0.01,0.01,1,0.01,0.01,0.01,0.01,1,0.01,0.01,0.01]
sig=[0.001,0.001,0.01,0.01,1,0.001,1,0.001,0.001]
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
np.savetxt(path+'B1422_flatchain_'+str(case_idx)+'.txt',sampler.flatchain)

