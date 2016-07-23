import numpy as np
import matplotlib.pyplot as plt

import DistanceTool as distance
import Basic_class as bClass

import pandas as pd

## read-in realization files & calculate f_sub w/i annulus Re(+0.05",-0.05")

## M_sub = pi*r_t*b_sub*Sigma_c

lensPath = '../../models/lens_info/'
lens_name = 'B1422'
gravlensPath = '../../data/sub_gravlens/'
realPath = gravlensPath+lens_name+'_realization_2/'

num = 1

filePath = realPath+'realization_'+str(num)+'.dat'

b_sub = np.loadtxt(filePath,dtype='float',unpack=True,usecols=[1])
x = np.loadtxt(filePath,dtype='float',unpack=True,usecols=[2])
y = np.loadtxt(filePath,dtype='float',unpack=True,usecols=[3])
rt = np.loadtxt(filePath,dtype='float',unpack=True,usecols=[9])

# get rid of macro model
b_sub,x,y,rt = b_sub[1:],x[1:],y[1:],rt[1:]


## Lens & Cosmology class
delta = 0.05 # arcsec

# Lens setup file
lens_setup = np.loadtxt(lensPath+lens_name+'_setup.dat')

Lens = bClass.Lens(lens_setup)

## cosmology parameter
cospara = bClass.Cosmology()

# critical surface density
sigma_c = Lens.critical_density()*cospara.h # M_sun/Mpc^2
sigma_c = sigma_c/(distance.mpc2arcs(cospara,1.,Lens.zl))**2 # M_sub/arcsec^2

## ------annulus mock ------ ##

r_out = Lens.b+delta
r_in = Lens.b-delta

r = np.sqrt(x**2+y**2)

# ----
cri_out = r < r_out

b_sub,x,y,rt = b_sub[cri_out],x[cri_out],y[cri_out],rt[cri_out]
r = r[cri_out]

# ---
cri_in = r > r_in

b_sub,x,y,rt = b_sub[cri_in],x[cri_in],y[cri_out],rt[cri_out]


## -----------

# calculate mass w/i rt


