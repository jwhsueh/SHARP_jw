import gravlens_tool as gt
import numpy as np

n_real = 10
n_start = 0

'''
## obs flux ratio
flux_ratio = np.array([1.062,	0.551,	0.024])
flux_err = np.array([0.009,	0.007,	0.006])
'''

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
optfile = 'B1422_real0_opt.input'

#####

## run for one set(360)

real_file = open(path+'/B1422_realization_nfw/flux_result_test.txt','w')
model_file = open(path+'/B1422_realization_nfw/smooth_result.txt','w')

for i in range(n_real):

	print '######'
	print 'realization '+str(i+n_start)
	print '######'

	table = np.loadtxt('../../data/sub_gravlens/B1422_realization_nfw/real'+str(i+n_start)+'.txt')
	x_list,y_list,ks_list,rs_list = table[:,1],table[:,2],table[:,3],table[:,4]


	class micro_mod:
		ks_i = ks_list
		rs_i = rs_list
		xi = x_list
		yi = y_list


	gt.create_opt(macro_mod,micro_mod,path,optfile)
	outline = gt.run_opt(path,optfile)

	outline = outline.split('\n')
	#print outline[21]
	#print outline[-18]
	macro_line = outline[21]
	src_line = outline[-18]
	chiline = outline[-17]
	#print chiline.split()
	macro_line = macro_line.split()[1:]
	src_line = src_line.split()[2:4]
	print macro_line,src_line
	chi2_pos, chi2_f = float(chiline.split()[2]),float(chiline.split()[3])
	print chi2_pos,chi2_f
	real_file.write(str(chi2_pos+chi2_f)+'\n')

	imgline = outline[-9:-5]

	#obs_x,obs_y,obs_f = np.empty(0),np.empty(0),np.empty(0)
	mod_x,mod_y,mod_f = np.empty(0),np.empty(0),np.empty(0)
	for line in imgline:
		elements = line.split()
		#obs_x = np.append(obs_x,float(elements[0]))
		#obs_y = np.append(obs_y,float(elements[1]))
		#obs_f = np.append(obs_f,float(elements[3]))

		mod_x = np.append(mod_x,float(elements[9]))
		mod_y = np.append(mod_y,float(elements[10]))
		mod_f = np.append(mod_f,float(elements[11]))


	## --- save mod_f and chi2
	print mod_f
	#real_file.write(str(mod_x)[1:-1]+'\t'+str(mod_y)[1:-1]+'\t'+str(mod_f)[1:-1]+'\t'+str(chi2_pos)+'\t'+str(chi2_f)+'\n')
	model_file.write(str(macro_line)[1:-2]+'\t'+str(src_line)[1:-1]+'\n')
	'''
	## check if img position has huge shift
	if np.sum((obs_x-mod_x)**2+(obs_y-mod_y)**2) > (3e-3)**2:
		print 'Discard this realization due to huge position shift'
	else:
		print 'Go to chi2 calculation'

		fr_mod = np.array([mod_f[1]/mod_f[0],mod_f[2]/mod_f[0],mod_f[3]/mod_f[0]])
		print fr_mod

		chi2 = (fr_mod-flux_ratio)**2/flux_err**2
		print chi2
		chi2 = np.sum(chi2)
		print chi2/3

		real_file.write(str(fr_mod)[1:-1]+'\t'+str(chi2/3)+'\n')

	'''
