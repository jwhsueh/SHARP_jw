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
		micro_line = gravlens_pjaffe(micro_mod.bi[i],micro_mod.xi[i],micro_mod.yi[i],micro_mod.rt_i[i])
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

	return

def run_opt(path,optfile):

	os.system('./lensmodel '+path+optfile)

	return

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

def create_findimg(micro_mod,path,output):
	findimg_file = open(path+output,'w')

	create_bestSub(path)
	opt_result = open(path+'bestSub.dat','r')

	## need to modify
	## write number of lens and srcs [startup]

	# get another function to create pure best substructure result

	lines = opt_result.readlines()[:-1]
	n_lens = len(lines)
	findimg_file.write('startup '+str(n_lens)+' 1\n')

	#opt_result.readline() # get rid of first line
	#opt_lines = opt_result.readlines()

	# write in substructures

	for Aline in lines:
		#Aline = opt_result.readline()
		findimg_file.write(Aline)
		#print Aline
		#print '##'+str(i)

	opt_result.close()

	## --end of gravlens startup--

	## --write a bunch of zeros--

	zeros = '0 0 0 0 0 0 0 0 0 0\n'
	for i in range(n_lens):
		findimg_file.write(zeros)

	findimg_file.write('\n')

	## write findimg command
	# get srcs position
	opt_result = open('best.dat','r')
	raw_line = opt_result.readlines() # get all the rest lines
	#print raw_line
	raw_line = raw_line[-11].split()
	#print raw_line

	findimg_file.write('findimg '+raw_line[1]+' '+raw_line[2]+'\n')

	opt_result.close()
	findimg_file.close()

	return

def gravlens_SIE(model):

	lenspara = [model.b,model.xc,model.yc,model.q,model.PA,model.gamma1,model.gamma2]
	lenspara = str(lenspara)[1:-1].replace(',','')
	Aline = 'alpha '+lenspara+' 0.0 0.0 1.0\n'

	return Aline

def gravlens_pjaffe(b,x,y,rt):

	lenspara = [b,x,y,0.0,0.0,0.0,0.0,0.0,rt,0.0]
	lenspara = str(lenspara)[1:-1].replace(',','')
	Aline = 'pjaffe '+lenspara[1:-1]+'\n'

	return Aline

def run_findimg(path,findimg_file):
	findimg_out = commands.getstatusoutput('./lensmodel '+path+findimg_file)

	findimg_out=findimg_out[1].split('\n')

	#print findimg_out

	result = open(path+'findimg.out','w')

	## checking img #

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
		print '* Not a quad-system'
		flag = False
		result.write('# Not a quad-system')

	result.close()

	return flag

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

