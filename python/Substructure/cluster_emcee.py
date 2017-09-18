import numpy as np
#import matplotlib.pyplot as plt
import gravlens_tool as gt
import emcee
#import sys

### --- findimg process

def do_findimg(paras,idx_sub,micro_mod,path_sub,path_macro,lens,idx_a,zlens,zsrc,burn_flag):
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
        zl = zlens
        zs = zsrc

    obs_file = path_macro+lens+'_obs.txt'
    table = np.loadtxt(obs_file)
    obs_x, obs_y, obs_f = table[:,0],table[:,1],table[:,2]
    err_x,err_y,err_f =  table[:,3],table[:,3],table[:,4]

    imgfile = 'real'+str(idx_sub)+'_findimg.input'
    gt.create_findimg(step_mod,micro_mod,path_sub,imgfile)
    qflag = gt.run_findimg(path_sub,imgfile,idx_sub) # flag tells you if it's quad

    if qflag==True:
        mod_x,mod_y,mod_f = gt.get_imgresult(path_sub,idx_sub)
        mod_x,mod_y,mod_f = gt.findimg_sort(mod_x,mod_y,mod_f)
        mod_f = np.abs(mod_f)
        mod_f = mod_f/mod_f[idx_a]
        print mod_f,obs_f
        print mod_x,obs_x
            
        chi2=0
        for i in range(4):
            chi2=chi2+(mod_x[i]-obs_x[i])**2/err_x[i]**2+(mod_y[i]-obs_y[i])**2/err_y[i]**2+(mod_f[i]-obs_f[i])**2/err_f[i]**2

        chi2 = chi2/9.0

        print chi2

        if (burn_flag!=0):
            #print mod_x,mod_y,mod_f
            chi2_outfile.write(str(chi2)+'\n')
            param_outfile.write(str(paras)[1:-1]+'\n')
    else:
        chi2=100000

    #print 'chi2 = '+str(chi2)

    return chi2

### --- log prob for emcee
def lnprob(paras,idx_sub,micro_mod,path_sub,path_macro,lens,idx_a,zlens,zsrc,burn_flag):
    #chi2 = do_findimg(paras,idx_sub,micro_mod,path_sub,path_macro,lens,idx_a,zlens,zsrc,burn_flag)

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
        zl = zlens
        zs = zsrc

    #obs_file = path_macro+lens+'_obs.txt'
    #table = np.loadtxt(obs_file)

    obs_x = np.array([-3.338751e-01,0,3.892583e-01,9.506678e-01])
    obs_y = np.array([-7.477099e-01,0, 3.199784e-01,-8.021654e-01])
    obs_f = np.array([0.551,1.062,1.0,0.024])
    err_x = np.array([0.00005,0.00005,0.00005,0.00005])
    err_y = np.array([0.00005,0.00005,0.00005,0.00005])
    err_f = np.array([0.007,0.009,1e-5,0.006])
    #obs_x, obs_y, obs_f = table[:,0],table[:,1],table[:,2]
    #err_x,err_y,err_f =  table[:,3],table[:,3],table[:,4]

    imgfile = 'real'+str(idx_sub)+'_findimg.input'
    gt.create_findimg(step_mod,micro_mod,path_sub,imgfile)
    qflag = gt.run_findimg(path_sub,imgfile,idx_sub) # flag tells you if it's quad

    if qflag==True:
        mod_x,mod_y,mod_f = gt.get_imgresult(path_sub,idx_sub)
        mod_x,mod_y,mod_f = gt.findimg_sort(mod_x,mod_y,mod_f)
        mod_f = np.abs(mod_f)
        mod_f = mod_f/mod_f[idx_a]
        print mod_f,obs_f
        print mod_x,obs_x
            
        chi2=0
        for i in range(4):
            chi2=chi2+(mod_x[i]-obs_x[i])**2/err_x[i]**2+(mod_y[i]-obs_y[i])**2/err_y[i]**2+(mod_f[i]-obs_f[i])**2/err_f[i]**2

        chi2 = chi2/9.0

        print chi2

        if (burn_flag!=0):
            #print mod_x,mod_y,mod_f
            chi2_outfile.write(str(chi2)+'\n')
            param_outfile.write(str(paras)[1:-1]+'\n')
    else:
        chi2=100000

    return -0.5*np.absolute(chi2)


### --- uniform random starting point

def start_point(mean,sig,w,nw):
# generate random start point
    
    ndim=len(mean)
    
    pos=np.zeros((nw,ndim))
    
    for i in range(nw):
        c=np.random.rand(ndim)
        for j in range(ndim):
            pos[i,j]=mean[j]-w*sig[j]+2*w*sig[j]*c[j]

    return pos


### -----   main function ------ ##

def emcee_sub(idx_sub,path_sub,path_macro,lens,zlens,zsrc,idx_a):

    ## read in sub realization
    table = np.loadtxt(path_sub+'real'+str(idx_sub)+'.txt')
    x_list,y_list,ks_list,rs_list = table[:,1],table[:,2],table[:,3],table[:,4]
    
    class micro_mod:
        ks_i = ks_list
        rs_i = rs_list
        xi = x_list
        yi = y_list
    '''
    class micro_mod:
        ks_i = np.zeros(1)
        rs_i = np.zeros(1)
        xi = np.zeros(1)
        yi = np.zeros(1)
    '''
    ## read in macro model starting point
    table = np.loadtxt(path_macro+lens+'_mean.txt')
    mean = table[0,:]

    ## define all output files
    #findimg_outfile = open(path_sub+'mod'+str(idx_sub)+'.txt','w')
    chi2_outfile = open(path_sub+'chi2_'+str(idx_sub)+'.txt','w')
    param_outfile = open(path_sub+'param_'+str(idx_sub)+'.txt','w')   

    ## --- emcee run settings

    nwalker=50 #number of chains
    ndim=9   #number of parameters
    burn=200   #number of (burn in) collect step
    nstep=1000  #number of MCMC steps

    burn_flag = 0

    ## set up for start points

    # expdisk
    #mean=[macro_mod.b, macro_mod.xc,macro_mod.yc,macro_mod.e,macro_mod.PA,macro_mod.gamma1,macro_mod.gamma2,macro_mod.sx,macro_mod.sy]
    #sig=[0.01,0.01,0.01,0.01,1,0.01,0.01,0.01,0.01,1,0.01,0.01,0.01]
    sig=[0.001,0.001,0.001,0.01,1,0.001,1,0.001,0.001]
    #sig=np.zeros(13)

    w_factor=1.0 # range of random start point


    # random generate star points
    #p0=np.zeros((nwalker,ndim))

    p0=start_point(mean,sig,w_factor,nwalker)

    ## sampling

    sampler=emcee.EnsembleSampler(nwalker,ndim,lnprob,args =[idx_sub,micro_mod,path_sub,path_macro,lens,idx_a,zlens,zsrc,burn_flag])

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


    ## --- save the flat chain
    #np.savetxt(path+'emcee_flatchain_'+str(case_idx)+'.txt',sampler.flatchain)

