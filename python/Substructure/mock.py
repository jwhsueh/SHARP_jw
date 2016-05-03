import numpy as np
import NFWprofile as NFW
import DistanceTool.py as dist

class cospara:
	OM = 0.27
	h = 0.71

class lenspara:
	zl = 0.6
	zs = 2.0
	b = 0.5

halo = NFW.set_halopara(cospara,lenspara)

## f_sat
sig_c = dist.critical_density(cospara,lenspara)
f_sat = 0.01

sig_tot = f_sat*sig_c/2.0 # total convergence of substructures

## different number of substructure
num_1 = 100
num_2 = 10

