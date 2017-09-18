import numpy as np
import gravlens_tool as gt
import os
#import cluster_emcee as mcmc

## this function do gravlens ray-tracing

def gravlens_rt(idx_sub,path_sub,path_out,lens,zlens,zsrc,idx_a):

### --- create findimg file from one macro model+one sub realization

	## read in sub realization
	table = np.loadtxt(path_sub+'real'+str(idx_sub)+'.txt')
	x_list,y_list,ks_list,rs_list = table[:,1],table[:,2],table[:,3],table[:,4]

	rs_list = np.sqrt(x_list**2+y_list**2)
	mask = np.logical_and(rs_list>0.8*0.75,rs_list<1.2*0.75)
	x_list,y_list,ks_list,rs_list = x_list[mask],y_list[mask],ks_list[mask],rs_list[mask]
	#print len(x_list)
	#print ks_list[0]

	class micro_mod:
		ks_i = ks_list
		rs_i = rs_list
		xi = x_list
		yi = y_list

	## read in macro model emsemble
	table = np.loadtxt(path_out+lens+'_eta.txt')
	b_list,xc_list,yc_list = table[:,0],table[:,1],table[:,2]
	e_list,pa_list,g1_list,g2_list = table[:,3],table[:,4],table[:,5],table[:,6]
	sx_list,sy_list = table[:,7],table[:,8]
	w_list = table[:,9] # weighting for importance sampling

	#N_macro = len(b_list)
	N_macro = 1

	## --- ray-trace through N macro model on one sub realization

	like_file = open(path_sub+'likelihood_'+str(idx_sub)+'.txt','w')
	
	for i in range(N_macro):
		print 'macro model', i

		class macro_mod:
		    b =  b_list[i]
		    xc = xc_list[i]
		    yc = yc_list[i]
		    e =  e_list[i]
		    PA = pa_list[i]
		    gamma1 =g1_list[i]
		    gamma2 = g2_list[i] ## shear PA
		    sx =sx_list[i]# src position
		    sy = sy_list[i]
		    zl = zlens
		    zs = zsrc

		findimg_outfile = 'real'+str(idx_sub)+'_findimg.input'


		gt.create_findimg(macro_mod,micro_mod,path_sub,findimg_outfile)

		weight = w_list[i]

		## do findimg, calculate chi2

		qflag = gt.run_findimg(path_sub,findimg_outfile,idx_sub)

		# read in obs file [in the x-cord increasing ordering. Flux should be respect to maximum (so does err_f)]
		obs_file = path_out+lens+'_obs.txt'
		table = np.loadtxt(obs_file)
		obs_x, obs_y, obs_f = table[:,0],table[:,1],table[:,2]
		#obs_f = obs_f/np.max(obs_f)
		err_x,err_y,err_f =  table[:,3],table[:,3],table[:,4]

		if qflag==True:
			mod_x,mod_y,mod_f = gt.get_imgresult(path_sub,idx_sub)
			mod_x,mod_y,mod_f = gt.findimg_sort(mod_x,mod_y,mod_f)
			mod_f = np.abs(mod_f)
			mod_f = mod_f/mod_f[idx_a]

			#print mod_f,obs_f
			#print mod_x,obs_x

			chi2=0
			for i in range(4):
				chi2=chi2+(mod_x[i]-obs_x[i])**2/err_x[i]**2+(mod_y[i]-obs_y[i])**2/err_y[i]**2+(mod_f[i]-obs_f[i])**2/err_f[i]**2

			chi2 = chi2/9.0

		else:
			chi2=100000

	    ## calculate likelihood
	    
		print chi2
		likelihood = np.exp(-0.5*np.absolute(chi2))*weight
		#print likelihood
		
		## record likelihood
		if likelihood!=0.0:
			like_file = open(path_sub+'likelihood_'+str(idx_sub)+'.txt','a')
			like_file.write(str(likelihood)+'\n')
			like_file.close()
			print likelihood


	 ## --- end of macro model loop
	like_file.close()
	#os.system('rm '+path_sub+findimg_outfile)
	#os.system('rm '+path_sub+'findimg'+str(idx_sub)+'.out')

	return 


