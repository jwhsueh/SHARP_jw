import numpy as np
import gravlens_tool as gt
import matplotlib.pyplot as plt
import NFWprofile as NFW

lens='B1422'
lens_path = '/Users/jwhsueh/Documents/SHARP_jw/models/B1422/sub_test/'
lens_dat = np.loadtxt(lens_path+lens+'_smooth.dat')
imgfile = lens+'_sub_findimg.input'

img_x,img_y,img_f = lens_dat[:,0],lens_dat[:,1],lens_dat[:,2]
f_idx =2
img_f = img_f/img_f[f_idx]
Re = 0.75
best_mod = np.array([7.466427e-01, 7.645018e-01, -6.730645e-01,3.782100e-01, 5.556261e+01, 1.473046e-01, 5.249075e+01])
src = np.array([3.902597e-01, -4.156308e-01])
xc,yc = best_mod[1],best_mod[2]

zl_g = 0.6 
zs_g = 3.12
zg = np.array([zl_g,zs_g])
sig_c = 5.75E+10 #[M_Sun h^-1 arcsec^-2] here!!!
sig_c = sig_c/0.704 #[M_Sun arcsec^-2] 

sub_mass = 1e7
sub_c200 = NFW.c200(zl_g,sub_mass)
sub_r200 = NFW.r200(zl_g,sub_mass)
sub_rs = (sub_r200/sub_c200).value
sub_ks = sub_mass/(4*np.pi*sub_rs**2*sig_c*(np.log(1+sub_c200)-sub_c200/(1+sub_c200)))

#print sub_ks,sub_rs

lens_type = ['sie','nfw']

res = 1e-2
box_size = 2.0*Re

x = np.linspace(xc-box_size,xc+box_size,np.floor(box_size*2/res)+1)
y = np.linspace(yc-box_size,yc+box_size,np.floor(box_size*2/res)+1)
gx,gy = np.meshgrid(x,y)

np.savetxt(lens_path+lens+'_cord_x.txt',gx)
np.savetxt(lens_path+lens+'_cord_y.txt',gy)

## collect closest distance

r0 = (gx-img_x[0])**2+(gy-img_y[0])**2
r1 = (gx-img_x[1])**2+(gy-img_y[1])**2
r2 = (gx-img_x[2])**2+(gy-img_y[2])**2
r3 = (gx-img_x[3])**2+(gy-img_y[3])**2

r = np.array([r0,r1,r2,r3])
r = np.sqrt(r)
r_min = np.min(r,axis=0)
r_min = r_min.flatten()
#print r.shape, r_min.shape

#mask = r_min<0.3
sub_x,sub_y = gx.flatten(),gy.flatten()

fa_frac,aa_arc = np.zeros(len(sub_x)),np.zeros(len(sub_x))

## ray-tracing loop

for i in range(len(sub_x)):
	if r_min[i]<0.4:
		print i

		sub_mod = np.array([sub_ks,sub_x[i],sub_y[i],sub_rs])
		rt_mod = np.array([best_mod,sub_mod])

		gt.create_findimg_macro(rt_mod,lens_type,zg,src,lens_path,imgfile)
		qflag = gt.run_findimg(lens_path,lens,imgfile) # flag tells you if it's quad

		if qflag==True:
			mod_x,mod_y,mod_f = gt.get_imgresult(lens_path,lens)
			#print mod_x
			#mod_x,mod_y,mod_f = findimg_sort(mod_x,mod_y,mod_f)
			mod_f = np.abs(mod_f)
			mod_f = mod_f/mod_f[f_idx]

			## calculate flux anomaly fraction & astrometric anomaly (in arcsec)

			flux_def = np.abs(mod_f-img_f)/img_f
			fa_frac[i] = np.max(flux_def)
			print fa_frac[i]

			pos_def = np.sqrt(np.abs(mod_x-img_x)**2+np.abs(mod_y-img_y)**2)
			aa_arc[i] = np.max(pos_def)


fa_frac = fa_frac.reshape(len(x),len(x))
aa_arc = aa_arc.reshape(len(x),len(x))

np.savetxt(lens_path+lens+'_07_fa.txt',fa_frac)
np.savetxt(lens_path+lens+'_07_aa.txt',aa_arc)

plt.imshow(fa_frac)
plt.show()
