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

num = 5000

f_sub = np.zeros(num)

## Lens & Cosmology class
delta = 0.05 # arcsec

# Lens setup file
lens_setup = np.loadtxt(lensPath+lens_name+'_setup.dat')

Lens = bClass.Lens(lens_setup)

r_out = Lens.b+delta
r_in = Lens.b-delta


## bring in the loop

for num_i in range(num):

	print num_i

	filePath = realPath+'realization_'+str(num_i+1)+'.dat'

	b_sub = np.loadtxt(filePath,dtype='float',unpack=True,usecols=[1])
	x = np.loadtxt(filePath,dtype='float',unpack=True,usecols=[2])
	y = np.loadtxt(filePath,dtype='float',unpack=True,usecols=[3])
	rt = np.loadtxt(filePath,dtype='float',unpack=True,usecols=[9])

# get rid of macro model
	b_sub,x,y,rt = b_sub[1:],x[1:],y[1:],rt[1:]

## ------annulus mock ------ ##

	area = np.pi*(r_out**2 - r_in**2)

## create a pandas table

	subs = pd.DataFrame()
	subs['b'],subs['x'],subs['y'],subs['r_t'] = b_sub,x,y,rt
	subs['r'] = np.sqrt(subs.x**2+subs.y**2)


	subs = subs.loc[(subs.r<r_out)&(subs.r>r_in),:]

## -----------

	# calculate mass w/i rt (in unit of Sigma_c)

	subs['m_rt'] = np.pi*subs.r_t*subs.b*0.586

	# substructure sigma

	subs['sigma'] = subs.m_rt/area

	f_sub[num_i] = np.sum(subs.sigma)*2


## plot & output

table = pd.DataFrame()
table['f_sub'] = f_sub

table.to_csv(realPath+'f_sub.dat',index = False)

table.hist(column = 'f_sub',bins = 20)
plt.title('B1422, f_sub w/i 0.1" annulus')
plt.xlabel('f_sub')
plt.ylabel('counts')
plt.savefig(realPath+'f_sub.png')



