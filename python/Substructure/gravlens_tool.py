import os
import numpy as np
import commands

def create_opt(macro_mod,micro_mod,path,output):
	opt_file = open(path+output,'w')
	opt_header = open(path+'opt.head','r')

	## write head file
	head_lines = opt_header.readlines()

	for Aline in head_lines:
		opt_file.write(Aline)

	## write number of lens and srcs [startup]

	n_lens = 1+len(micro_mod.xi)

	opt_file.write('startup '+str(n_lens)+' 1\n')

	## macro model
	macro_line = gravlens_SIE(macro_mod)
	opt_file.write(macro_line)

	## micro model
	for i in range(len(micro_mod.xi)):
		#micro_line = gravlens_pjaffe(micro_mod.bi[i],micro_mod.xi[i],micro_mod.yi[i],micro_mod.rt_i[i])
		micro_line = gravlens_nfw(micro_mod.ks_i[i],micro_mod.xi[i],micro_mod.yi[i],micro_mod.rs_i[i])
		opt_file.write(micro_line)

	## --end of gravlens startup--

	## write opt SIE flags

	opt_file.write('\n')
	opt_file.write('1 1 1 1 1 1 1 0 0 0\n')

	## --write a bunch of zeros [pjaffe]--

	
	zeros = '0 0 0 0 0 0 0 0 0 0\n'
	for i in range(n_lens-1):
		opt_file.write(zeros)

	opt_file.write('\n')

	## write opt command

	opt_file.write('optimize \n')

	opt_file.close()



def create_sieopt(macro_mod,path,datafile,filename):

	opt_file = open(path+filename,'w')

	opt_file.write('set omitcore=1.0e-6\n')
	opt_file.write('data '+path+datafile+'\n')
	opt_file.write('set checkparity = 0\nset gridflag= 0\nstartup 1 1\n')

	## macro model
	macro_line = gravlens_SIE(macro_mod)
	opt_file.write(macro_line)

	## --end of gravlens startup--

	## write opt SIE flags

	opt_file.write('\n')
	opt_file.write('1 1 1 1 1 1 1 0 0 0\n')

	## write opt command

	opt_file.write('set chimode= 0\noptimize \nset chimode= 1\noptimize \noptimize \n')

	opt_file.close()

def gravlens_opt(path,opt_filename):

## run opt 
		
	opt_out = commands.getstatusoutput('./lensmodel '+path+opt_filename)
	opt_out = opt_out[1].split('\n')
	#print opt_out
	print opt_out[-5]

	## read best.dat
	# img
	bestfile = open('best.dat','r')
	best_line = []
	for line in bestfile.readlines():
		best_line.append(line)
	##

	lens_model = best_line[1]
	lens_para = lens_model.split()
	lens_para = np.array(lens_para[1:8]).astype(float)
	print lens_para
	#shear = float(lens_para[6])
	#len_x,len_y = float(lens_para[2]),float(lens_para[3])
	#len_cen = np.array([len_x,len_y])
	#print 'shear= '+str(shear)

	# chi2

	chi2 = best_line[4].split()
	chi2 =float(chi2[2])

	print 'chi-square='+str(chi2)

	# src
	src = best_line[7].split()
	src = np.array([float(src[1]),float(src[2])])
	print src

	return best_line,lens_para,src,chi2

	

def run_opt(path,optfile):

	outline = commands.getstatusoutput('./lensmodel '+path+optfile)

	return outline[1]

def create_bestSub(path):

	best_out = open('best.dat','r')

	lines = best_out.readlines()

	lines = lines[1:-16]
	#print lines

	best_sub = open(path+'bestSub.dat','w')

	for Aline in lines:
		best_sub.write(Aline)

	best_out.close()
	best_sub.close()

	return

def create_findimg_micro(macro_mod,micro_mod,path,output):
	findimg_file = open(path+output,'w')

	## write zl, zs
	findimg_file.write('set zlens='+str(macro_mod.zl)+'\n')
	findimg_file.write('set zsrc='+str(macro_mod.zs)+'\n')

	len2 = len(micro_mod.xi)
	n_lens = 1+len2
	findimg_file.write('startup '+str(n_lens)+' 1\n')

	## --- set up lens
	# macro model
	## here
	model_line = np.array([macro_mod.b,macro_mod.xc,macro_mod.yc,macro_mod.e,macro_mod.PA,macro_mod.gamma1,macro_mod.gamma2])
	macro_line = gravlens_SIE(model_line)
	findimg_file.write(macro_line)

	## micro model
	for i in range(len(micro_mod.xi)):
		#micro_line = gravlens_pjaffe(micro_mod.bi[i],micro_mod.xi[i],micro_mod.yi[i],micro_mod.rt_i[i])
		micro_line = gravlens_nfw(micro_mod.ks_i[i],micro_mod.xi[i],micro_mod.yi[i],micro_mod.rs_i[i])
		findimg_file.write(micro_line)

	## --end of gravlens startup--

	## --write a bunch of zeros--

	zeros = '0 0 0 0 0 0 0 0 0 0\n'
	for i in range(n_lens):
		findimg_file.write(zeros)

	findimg_file.write('\n')

	## write findimg command
	# get srcs position

	findimg_file.write('findimg '+str(macro_mod.sx)+' '+str(macro_mod.sy)+'\n')

	findimg_file.close()

	return 

def create_findimg_macro(macro_mod,profile,z,src,path,output):

	findimg_file = open(path+output,'w')
	#print "create findimg file"

	zl,zs = z[0],z[1]

	## write zl, zs
	findimg_file.write('set zlens='+str(zl)+'\n')
	findimg_file.write('set zsrc='+str(zs)+'\n')

	n_lens = macro_mod.shape[0]

	findimg_file.write('startup '+str(n_lens)+' 1\n')
	#print 'here'

	## --- set up lens
	# macro model
	#print macro_mod
	for i in range(n_lens): ## here
		if profile[i] == 'sie':
			macro_line = gravlens_SIE(macro_mod[i])
			findimg_file.write(macro_line)
		if profile[i] == 'nfw':
			macro_line = gravlens_nfw(macro_mod[i])
			findimg_file.write(macro_line)

	## --write a bunch of zeros--

	zeros = '0 0 0 0 0 0 0 0 0 0\n'
	for i in range(n_lens):
		findimg_file.write(zeros)

	findimg_file.write('\n')

	## write findimg command
	# get srcs position

	findimg_file.write('findimg '+str(src[0])+' '+str(src[1])+'\n')

	findimg_file.close()

	return 

def create_caustic(macro_model,lens_z,infile,outfile,path):

	gravlens_file = open(path+infile,'w')

	zl,zs = lens_z[0],lens_z[1]

	## write zl, zs
	gravlens_file.write('set zlens='+str(zl)+'\n')
	gravlens_file.write('set zsrc='+str(zs)+'\n')

	gravlens_file.write('startup 1 1\n')
	#print 'here'

	## --- set up lens
	# macro model

	macro_line = gravlens_SIE(macro_model)
	gravlens_file.write(macro_line)

	zeros = '0 0 0 0 0 0 0 0 0 0\n'
	gravlens_file.write(zeros)

	gravlens_file.write('\n')

	## --- plot caustic
	gravlens_file.write('plotcrit '+path+outfile)

#def run_gravlens_file(filename,path):
	
def quad_src(n_src,lens_cen,path,filename):

	tab = np.loadtxt(path+filename)
	ux,uy = np.append(tab[:,2],tab[:,6]),np.append(tab[:,3],tab[:,7])

	# dist to lens center
	u_dist = np.sqrt((ux-lens_cen[0])**2+(uy-lens_cen[1])**2)
	u_mid = np.average(u_dist)

	# 1st & 3rd quad
	u_sm = u_dist[u_dist<=u_mid]
	u_la = u_dist[u_dist>u_mid]
	u1,u3 = np.average(u_sm),np.average(u_la)

	# grouping -> tangential caustic
	mask = np.abs(u_dist-u1) <= np.abs(u_dist-u3)
	ux,uy,ur = ux[mask],uy[mask],u_dist[mask]

	x0,x1,y0,y1 = np.min(ux),np.max(ux),np.min(uy),np.max(uy)
	r_min = np.min(u_dist)

	sx0,sy0 = 0.0,0.0
	src_x,src_y = [],[] # collect src pos
	i = 0
	while i<n_src:
		## randomly sample
		xi = np.random.random_sample()*(x1-x0)+x0
		yi = np.random.random_sample()*(y1-y0)+y0

		# distance to caustic
		ri = np.sqrt((xi-ux)**2+(yi-uy)**2)
		rci = np.sqrt((xi-lens_cen[0])**2+(yi-lens_cen[1])**2) # dist to center
		r_dist = np.min(ri) # shortest dist to caustic
		idx = np.argmin(ri)
		rc = ur[idx] # dist to center of that caus pt

		if np.logical_and(r_dist<0.25*r_min,r_dist>0.01*r_min):
			if (rci<rc):
				if np.logical_and(xi!=sx0,yi!=sy0):
					src_x.append(xi)
					src_y.append(yi)

					sx0,sy0=xi,yi
					i = i+1

	return np.array(src_x),np.array(src_y)

def create_mockdat(img_x,img_y,img_f,err_x,path,filename):
	datfile = open(path+filename,'w')

	datfile.write('1\n')
	datfile.write('0 0 100000\n')
	datfile.write('0.3 100000\n')
	datfile.write('99 100000\n')
	datfile.write('0.2 100000\n\n')
	datfile.write('1\n4\n')

	datfile.write(str(img_x[0])+'\t'+str(img_y[0])+'\t'+str(img_f[0])+'\t'+str(err_x)+' 1000  0.0 0.0 \n')
	datfile.write(str(img_x[1])+'\t'+str(img_y[1])+'\t'+str(img_f[1])+'\t'+str(err_x)+' 1000  0.0 0.0 \n')
	datfile.write(str(img_x[2])+'\t'+str(img_y[2])+'\t'+str(img_f[2])+'\t'+str(err_x)+' 1000  0.0 0.0 \n')
	datfile.write(str(img_x[3])+'\t'+str(img_y[3])+'\t'+str(img_f[3])+'\t'+str(err_x)+' 1000  0.0 0.0 \n')

	datfile.close()

def gravlens_SIE(model):

	lenspara = [model[0],model[1],model[2],model[3],model[4],model[5],model[6]] 
	# model.b,model.xc,model.yc,model.e,model.PA,model.gamma1,model.gamma2
	lenspara = str(lenspara).replace(',','')
	Aline = 'alpha '+lenspara[1:-1]+' 0.0 0.0 1.0\n'

	return Aline

def gravlens_pjaffe(b,x,y,rt):

	lenspara = [b,x,y,0.0,0.0,0.0,0.0,0.0,rt,0.0]
	lenspara = str(lenspara).replace(',','')
	Aline = 'pjaffe '+lenspara[1:-1]+'\n'

	return Aline

def gravlens_nfw(model):

	ks,x,y,rs = model[0],model[1],model[2],model[3]

	lenspara = [ks,x,y,0.0,0.0,0.0,0.0,rs,0.0,0.0]
	lenspara = str(lenspara).replace(',','')
	Aline = 'nfw '+lenspara[1:-1]+'\n'

	return Aline

def run_findimg(path,lens,findimg_file):
	findimg_out = commands.getstatusoutput('./lensmodel '+path+findimg_file)

	findimg_out=findimg_out[1].split('\n')

	#print findimg_out

	result = open(path+lens+'_findimg.out','w')

	## checking img #
	#print findimg_out
	check_line = findimg_out[-6].split()
	# if it's four images
	if check_line[0] == '#':
		# save x y mag
		#x,y,mag = [],[],[]
		for i in [-2,-3,-4,-5]:
			Aline = findimg_out[i].split()
			#x.append(float(Aline[0]))
			#y.append(float(Aline[1]))
			#mag.append(float(Aline[2]))
			result.write(Aline[0]+' '+Aline[1]+' '+Aline[2]+'\n')

			flag = True

		#result = [x,y,mag] # write into a file


	else:
		#print '* Not a quad-system'
		flag = False
		result.write('# Not a quad-system')

	result.close()

	return flag

def get_imgresult(path,filename):
	table = np.loadtxt(path+filename+'_findimg.out')

	x,y,f = table[:,0],table[:,1],table[:,2]

	return x,y,f

def findimg_sort(mod_x,mod_y,mod_f):
    xsort = np.sort(mod_x)

    x_new,y_new,f_new = np.empty(4),np.empty(4),np.empty(4)
    for i in range(4):
        idx = list(mod_x).index(xsort[i])

        x_new[i],y_new[i],f_new[i] = mod_x[idx],mod_y[idx],mod_f[idx]

    return x_new,y_new,f_new

def create_lensclass(paras):
	class create_mod:
		b = paras[0]
    	xc = paras[1]
    	yc = paras[2]
    	e = paras[3]
    	PA = paras[4]
    	gamma1 = paras[5]
    	gamma2 = paras[6]
    	sx = paras[7]
    	sy = paras[8]
	print create_mod
	return create_mod

def assign_img(path,x_obs):
	img_file = path+'findimg.out'
	
	table = np.loadtxt(img_file)

	result = open(path+'findimg.out','w')

	for i in range(len(x_obs)):
		d = table[:,0] - x_obs[i] # difference in pos_x
		d = np.abs(d)
		d_min = np.min(d)
		idx = list(d).index(d_min)

		result.write(str(table[idx,:])[1:-1]+'\n')

	return

## ---- check if result is w/i 2-sigma----
def pos_check(path,obs_x,obs_y,obs_err):

	img_file = path+'findimg.out'
	
	table = np.loadtxt(img_file)
	#print obs_err
	mod_x,mod_y = table[:,0],table[:,1]
	dif_x = np.abs(obs_x - mod_x)
	dif_y = np.abs(obs_y - mod_y)
	obs_err = np.average(obs_err)

	print obs_x
	print mod_x
	print table[:,2]
	#print obs_err
	#print dif_x,dif_y

	result = True 

	#if np.average(dif_x) > 3.*obs_err or np.average(dif_y) > 3.*obs_err:
	#	result = False

	for element in dif_x:
		if element>2.*obs_err:
			result=False

	for element in dif_x:
		if element>2.*obs_err:
			result=False

	return result

## ----R_cusp, R_fold-----

def write_Rcusp(Rfile,path):
	mod_file = path+'findimg.out'

	flux = np.loadtxt(mod_file,dtype='float',unpack=True,usecols=[2])

	R_cusp = (flux[0]-flux[1]+flux[2])/np.abs(flux[0]+flux[1]+flux[2])

	Rfile.write(str(R_cusp)+'\n')

	return

def write_Rfold(Rfile,path):
	mod_file = path+'findimg.out'

	flux = np.loadtxt(mod_file,dtype='float',unpack=True,usecols=[2])

	R_fold = (flux[0]-flux[1])/np.abs(flux[0]+flux[1])

	Rfile.write(str(R_fold)+'\n')

	return

def write_chi2(chi_file,lens,path):
	mod_file = path+'findimg.out'
	x_mod = np.loadtxt(mod_file,dtype='float',unpack=True,usecols=[0])
	y_mod = np.loadtxt(mod_file,dtype='float',unpack=True,usecols=[1])
	f_mod = np.loadtxt(mod_file,dtype='float',unpack=True,usecols=[2])

	# flux ratio
	fr_mod = f_mod/f_mod[0]

	x_obs,y_obs,fr_obs = lens.img_x,lens.img_y,lens.img_fr
	pos_err,fr_err = lens.img_err,lens.img_frerr

	print x_obs,y_obs,fr_obs,pos_err,fr_err
	print x_mod,y_mod,fr_mod

	chi2 = 0.
	for i in range(x_obs.size):
		chi2 = chi2+(x_obs[i]-x_mod[i])**2/pos_err[i]**2+(y_obs[i]-y_mod[i])**2/pos_err[i]**2

		if i >0:
			chi2 = chi2 + (fr_obs[i]-fr_mod[i])**2/fr_err[i]**2


	chi_file.write(str(chi2)+'\n')

def write_realization(Rnum,path,Rpath):

	info = open(path+'bestSub.dat','r')
	real_file = open(Rpath+'realization_'+str(Rnum)+'.dat','w')

	lines = info.readlines()

	for Aline in lines:
		real_file.write(Aline)

	info.close()
	real_file.close()


### -------

def get_bestdat(filename):

	bestfile = open(filename,'r')
	lines = bestfile.readlines()
	#print lines[-6:]

	datafit = lines[-6:-2]
	img_x,img_y,img_f = [],[],[]
	for Aline in datafit:
		fitline =  Aline.split()
		img_x.append(float(fitline[-4]))
		img_y.append(float(fitline[-3]))
		img_f.append(float(fitline[-2]))


	img_x,img_y,img_f = np.array(img_x),np.array(img_y),np.array(img_f)
	return img_x,img_y,img_f
