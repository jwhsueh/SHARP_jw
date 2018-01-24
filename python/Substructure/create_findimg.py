import numpy as np
import gravlens_tool as gt
import sys

case_idx = int(sys.argv[1])
### --- create findimg file from one macro model+one sub realization

path = '../../models/mock1/gravlens_code/'

## read in sub realization
table = np.loadtxt('/Volumes/sting_1/subs/mock1/real_20/real'+str(case_idx)+'.txt')
x_list,y_list,ks_list,rs_list = table[:,1],table[:,2],table[:,0],table[:,3]

class macro_mod: ## B1422
    b =  2.086805e-01
    xc =2.104428e-01
    yc =  -1.067919e-01 
    e =  0.34
    PA = 2.426760e+01
    gamma1 =2.397602e-02 
    gamma2 = -4.430349e+01 ## shear PA
    sx =0.226   # src position
    sy = -0.1415
    zl = 0.6
    zs = 3.12

class micro_mod:
	ks_i = ks_list
	rs_i = rs_list
	xi = x_list+macro_mod.xc
	yi = y_list+macro_mod.yc

findimg_outfile = 'mock1_'+str(case_idx)+'_findimg.input'


gt.create_findimg(macro_mod,micro_mod,path,findimg_outfile)
