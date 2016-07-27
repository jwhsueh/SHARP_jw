""" This code generate realizations under fixed f_sub and call gravlens findimg """

import numpy as np
import matplotlib.pyplot as plt

import MassFunc as MF
import DistanceTool as distance
import Basic_class as bClass
import gravlens_tool as gTool

## initiate lens
lens_name = 'B1422'
lens_type = 'cusp'
lensPath = '../../models/lens_info/'
gravlensPath = '../../data/sub_gravlens/'
realPath = gravlensPath+lens_name+'_realization_2/'

# Lens setup file
lens_setup = np.loadtxt(lensPath+lens_name+'_setup.dat')
lens_obs = np.loadtxt(lensPath+lens_name+'_obs.dat')

#B1422 = bClass.Lens(z_lens= 0.536,z_src = 3.62,ellipticity = 0.4, Re = 0.7511, x_cen = 0.74,y_cen = -0.66)
Lens = bClass.Lens(lens_setup)
Lens.obsData(lens_obs)

## cosmology parameter
cospara = bClass.Cosmology()

## realization boundary
x_lim = np.array([Lens.xc-2.*Lens.b,Lens.xc+2.*Lens.b])
y_lim = np.array([Lens.yc-2.*Lens.b,Lens.yc+2.*Lens.b])


## substructure mass fraction
f_sub = 0.01

sigma_c = Lens.critical_density()*cospara.h # M_sun/Mpc^2
sigma_c = sigma_c/(distance.mpc2arcs(cospara,1.,Lens.zl))**2
sig_tot = f_sub*sigma_c/2.0 # M_sun/arcsec^2
print sig_tot


## 2D position + mass
## truncation radius (fixed at 1/60 -- F&K)

def uniform_2d(n_draw):
	x = np.random.uniform(x_lim[0],x_lim[1],n_draw)
	y = np.random.uniform(y_lim[0],y_lim[1],n_draw)

	dx = x - Lens.xc
	dy = y - Lens.yc

	r = np.sqrt(dx**2+dy**2)

	x = x[r<=2.*Lens.b]
	y = y[r<=2.*Lens.b]
	#r = r[r<=2.*Lens.b]

	return x,y

def get_mass(n_draw):
	draws = np.random.rand(n_draw)

	m_min,m_max = 1e6, 1e8 # lower & upper limits

	MF_x,MF_y = MF.inverse_cdf(m_min,m_max)

	m = np.interp(draws,MF_x,MF_y)

	return m

def m2b_jaffe(m_tot):
	# M_sub = pi*b*b_sub^3/2 * sigma_c 

	#b_mpc = distance.arcs2mpc(cospara,Lens.b,Lens.zl)

	bs_15 = m_tot/np.pi/Lens.b/sigma_c
	bs = np.power(bs_15,2./3.)
	#bs = distance.mpc2arcs(cospara,bs_mpc,Lens.zl)

	return bs

def jaffe_sig(m_sub): ## calculate sigma: m_sub, circle of r = 2*Re

	js = m_sub/(np.pi*(2.*Lens.b)**2)

	return js

## --- realization class (store realization info)

class realization_ob:
	def __init__(self,x,y,m,b,rt):
		self.xi,self.yi,self.mi,self.bi,self.rt_i = x,y,m,b,rt

## Function: create realization

def set_realization(): # need to add realization number

##---- Here we start the drawing ----##

	x_sub = []
	y_sub = []
	m_sub,b_sub,rt_sub = [],[],[]

	sig_idx = 0.

	unit_draw = 10 # 1 unit for the drawing

##----keep drawing substructure
	while True:
		
		print unit_draw
	# position

		xi,yi = uniform_2d(unit_draw)

		unit_sub = xi.size # number of substructure in this unit draw

		mi = get_mass(unit_sub)
		print mi

	#-----
	# total mass
	# mass to pseudo-Jaffe profile
		b_i = m2b_jaffe(mi)
		print b_i
		rt_i = np.sqrt(b_i*Lens.b)

	# sub sig calculation
		sig_sum = np.sum(mi)/np.pi/(2.*Lens.b)**2

	## set up counting approach [for loop]
		print sig_sum
		sig_idx = sig_idx+sig_sum

		if sig_idx < sig_tot:
			x_sub.extend(xi)
			y_sub.extend(yi)
			m_sub.extend(mi)
			b_sub.extend(b_i)
			rt_sub.extend(rt_i)

		else:
			sig_idx = sig_idx - sig_sum

			delta_sig = sig_tot - sig_idx

			xi,yi = xi[0],yi[0]
			mi = delta_sig*np.pi*(2.*Lens.b)**2
			b_i = m2b_jaffe(mi)
			rt_i = np.sqrt(b_i*Lens.b)

			x_sub.append(xi)
			y_sub.append(yi)
			m_sub.append(mi)
			b_sub.append(b_i)
			rt_sub.append(rt_i)

			break

	## ------end of drawing----			

	# create realization class object

	real_one = realization_ob(x_sub,y_sub,m_sub,b_sub,rt_sub)
		
	return real_one


## -------------------------------
valid_num = 5000
sucess = 0

cusp_file = open(gravlensPath+str(lens_name)+'_Rcusp_2.dat','a')
cusp_file.write('# f_sub='+str(f_sub)+', realization='+str(valid_num)+'\n')

fold_file = open(gravlensPath+str(lens_name)+'_Rfold_2.dat','a')
fold_file.write('# f_sub='+str(f_sub)+', realization='+str(valid_num)+'\n')

chi_file = open(gravlensPath+str(lens_name)+'_chi2_2.dat','a')
chi_file.write('# f_sub='+str(f_sub)+', realization='+str(valid_num)+'\n')


while True:
	re_mod = set_realization()
	print len(re_mod.xi)

	opt_name = 'opt_realization.input'
	findimg_name = 'findimg_realization.input'

	gTool.create_opt(Lens,re_mod,gravlensPath,opt_name)
	## run opt file
	gTool.run_opt(gravlensPath,opt_name)

	print '----Optimization done----'

	## write findimg file
	gTool.create_findimg(re_mod,gravlensPath,findimg_name)

	print '----create findimg file done----'

	flag = gTool.run_findimg(gravlensPath,findimg_name)

	print '----write findimg result----'

	

	if flag == False:
	# run above again
		print 'realization fail'

	else:
		gTool.assign_img(gravlensPath,Lens.img_x)
		check = gTool.pos_check(gravlensPath,Lens.img_x,Lens.img_y,Lens.img_err)

		if check == False:
			print '# realization fail'

		else: ##------created
			print '# realization created'
			sucess = sucess+1

			gTool.write_Rcusp(cusp_file,gravlensPath)
			gTool.write_Rfold(fold_file,gravlensPath)
			gTool.write_chi2(chi_file,Lens,gravlensPath)

			## save realization info
			gTool.write_realization(sucess,gravlensPath,realPath)







	if sucess>=valid_num:
		break


print 'final sucess = '
print sucess
