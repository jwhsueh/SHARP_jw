import numpy as np
import matplotlib.pyplot as plt
import gravlens_tool as gt
import emcee
import multiprocessing

class macro_mod: ## B1422
	b =  7.466419e-01
	xc =7.645039e-01
	yc =  -6.730659e-01
	e =  3.782159e-01
	PA = 5.556257e+01
	gamma1 =1.473030e-01
	gamma2 = 5.249074e+01 ## shear PA
	sx =3.902599e-01 # src position
	sy = -4.156310e-01

path = '../../data/sub_gravlens/'
imgfile = 'B1422_real0_findimg2.input' #findimg

table = np.loadtxt('../../data/sub_gravlens/B1422_realization1/real0.txt')
x_list,y_list,bsub_list,rt_list = table[:,1],table[:,2],table[:,3],table[:,4]


class micro_mod:
	bi = bsub_list
	rt_i = rt_list
	xi = x_list
	yi = y_list

obs_x = np.array([-3.338751e-01,1.629951e-06,3.892583e-01,9.506678e-01])
obs_y = np.array([-7.477099e-01,-8.062674e-06, 3.199784e-01,-8.021654e-01])
obs_f = np.array([3.936236e+00,8.322586e+00, 6.970661e+00,2.507707e-01])

err_x = np.array([0.001,0.001,0.001,0.001])
err_y = np.array([0.001,0.001,0.001,0.001])
err_f = obs_f*0.05

burn_flag = 0

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
    gt.create_findimg(step_mod,micro_mod,path,imgfile)
    qflag = gt.run_findimg(path,imgfile) # flag tells you if it's quad

    if qflag==True:
        mod_x,mod_y,mod_f = gt.get_imgresult(path)
        mod_x,mod_y,mod_f = findimg_sort(mod_x,mod_y,mod_f)
        
        if (burn_flag!=0):
            print mod_x,mod_y,mod_f
            
        chi2=0
        for i in range(4):
            chi2=chi2+(mod_x[i]-obs_x[i])**2/err_x[i]**2+(mod_y[i]-obs_y[i])**2/err_y[i]**2+(mod_f[i]-obs_f[i])**2/err_f[i]**2

        chi2 = chi2/9.0
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
burn=20   #number of (burn in) collect step
nstep=1000  #number of MCMC steps

## set up for start points

# expdisk
mean=[macro_mod.b, macro_mod.xc,macro_mod.yc,macro_mod.e,macro_mod.PA,macro_mod.gamma1,macro_mod.gamma2,macro_mod.sx,macro_mod.sy]
#sig=[0.01,0.01,0.01,0.01,1,0.01,0.01,0.01,0.01,1,0.01,0.01,0.01]
sig=[0.002,0.001,0.001,0.01,1,0.01,1,0.001,0.001]
#sig=np.zeros(13)

w_factor=1.0 # range of random start point


# random generate star points
#p0=np.zeros((nwalker,ndim))

p0=start_point(mean,sig,w_factor,nwalker)

## sampling

sampler=emcee.EnsembleSampler(nwalker,ndim,lnprob,threads=15)

## burn-in stage
print '### burn-in stage started ###'
print '\n'

pos, prob, state = sampler.run_mcmc(p0,burn)

burn_flag = 1
print '### burn-in stage finished ###'
print '\n'
print '### star MCMC run ###'

## MCMC run
sampler.reset()
sampler.run_mcmc(pos,nstep)
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


