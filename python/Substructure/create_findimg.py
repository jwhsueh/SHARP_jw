import numpy as np
import gravlens_tool as gt
import sys

case_idx = int(sys.argv[1])
### --- create findimg file from one macro model+one sub realization

path = '../../data/sub_gravlens/'

## read in sub realization
table = np.loadtxt('../../data/sub_gravlens/B1422_realization1/real'+str(case_idx)+'.txt')
x_list,y_list,bsub_list,rt_list = table[:,1],table[:,2],table[:,3],table[:,4]

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

class micro_mod:
	bi = bsub_list
	rt_i = rt_list
	xi = x_list
	yi = y_list

findimg_outfile = 'mock'+str(case_idx)+'_findimg.input'


gt.create_findimg(macro_mod,micro_mod,path,findimg_outfile)
