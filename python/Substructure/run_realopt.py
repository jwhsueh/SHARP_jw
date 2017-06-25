import gravlens_tool as gt
import numpy as np

## obs flux ratio
flux_ratio = np.array([1.062,	0.551,	0.024])
flux_err = np.array([0.009,	0.007,	0.006])

class macro_mod: ## B0128
	b =  7.526814e-01 
	xc = 7.486039e-01
	yc =  -6.622528e-01
	e = 3.323503e-01
	PA = 5.594659e+01
	gamma1 =1.595633e-01
	gamma2 = 5.256553e+01 ## shear PA

path = '../../data/sub_gravlens/'
optfile = 'B1422_real0_opt.input'

#####

## run for one set(360)

real_file = open(path+'/B1422_realization/flux_result.txt','w')

for i in range(360):

	print '######'
	print 'realization '+str(i)
	print '######'

	table = np.loadtxt('../../data/sub_gravlens/B1422_realization/real'+str(i)+'.txt')
	x_list,y_list,bsub_list,rt_list = table[:,1],table[:,2],table[:,3],table[:,4]


	class micro_mod:
		bi = bsub_list
		rt_i = rt_list
		xi = x_list
		yi = y_list


	gt.create_opt(macro_mod,micro_mod,path,optfile)
	outline = gt.run_opt(path,optfile)

	outline = outline.split('\n')
	imgline = outline[-9:-5]

	obs_x,obs_y,obs_f = np.empty(0),np.empty(0),np.empty(0)
	mod_x,mod_y,mod_f = np.empty(0),np.empty(0),np.empty(0)
	for line in imgline:
		elements = line.split()
		obs_x = np.append(obs_x,float(elements[0]))
		obs_y = np.append(obs_y,float(elements[1]))
		obs_f = np.append(obs_f,float(elements[3]))

		mod_x = np.append(mod_x,float(elements[9]))
		mod_y = np.append(mod_y,float(elements[10]))
		mod_f = np.append(mod_f,float(elements[11]))
		#print elements

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


