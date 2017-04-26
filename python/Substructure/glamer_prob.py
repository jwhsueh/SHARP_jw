import numpy as np

path='/Volumes/sting_1/data/'

list_file = path+'snap99_fa_Rcusp.txt'

file_list = np.genfromtxt(list_file,dtype='str')

Rfold,Rcusp,phi0,phi1 = np.empty(1),np.empty(1),np.empty(1),np.empty(1)
for file_one in file_list:
	table = np.loadtxt(path+'Rcusp_f/'+file_one+'_64_Rcusp_ga.txt')
	Rfold,Rcusp,phi0,phi1 = np.append(Rfold,table[:,0]),np.append(Rcusp,table[:,1]),np.append(phi0,table[:,2]),np.append(phi1,table[:,3])

mask = np.abs(Rfold)<0.7
Rfold,Rcusp,phi0,phi1=Rfold[mask],Rcusp[mask],phi0[mask],phi1[mask]
mask = np.abs(Rcusp)<0.7
Rfold,Rcusp,phi0,phi1=Rfold[mask],Rcusp[mask],phi0[mask],phi1[mask]

## ---- real lens data
# elp
tab_e = np.loadtxt('../../data/flux_ratio.txt')
Rf_e,Rc_e,p0_e,p1_e,Rfe_e,Rce_e = tab_e[:,3],tab_e[:,1],tab_e[:,0],tab_e[:,2],tab_e[:,4],tab_e[:,5]

tab_d = np.loadtxt('../../data/flux_ratio_disk.txt')
Rf_d,Rc_d,p0_d,p1_d,Rfe_d,Rce_d = tab_d[:,3],tab_d[:,1],tab_d[:,0],tab_d[:,2],tab_d[:,4],tab_d[:,5]

p0,p1,Rfe,Rce=p0_e,p1_e,Rfe_e,Rce_e
Rf,Rc=np.abs(Rf_e),np.abs(Rc_e)

#p0,p1,Rfe,Rce=p0_d,p1_d,Rfe_d,Rce_d
#Rf,Rc=np.abs(Rf_d),np.abs(Rc_d)

## ---- probability contour

## cusp
print "#### cusp ####"

for i in range(len(p0)):
	pcen,Rc_real = p0[i],Rc[i]
	Rc_p,Rc_m = Rc[i]+Rce[i],Rc[i]-Rce[i]
	print pcen,Rc_real
	mask = np.logical_and(phi0>=pcen-5,phi0<pcen+5)
	
	Rcusp_b = np.abs(Rcusp[mask])
	print Rcusp_b.size
	sig = np.sqrt(Rcusp_b.size)

	Rcusp_b = np.sort(Rcusp_b)
	delta = np.abs(Rcusp_b-Rc_real)
	idx = np.argmin(delta)
	prob = 1.0-float(idx)/Rcusp_b.size

	delta = np.abs(Rcusp_b-Rc_p)
	idx = np.argmin(delta)
	prob_p = 1.0-float(idx)/Rcusp_b.size

	delta = np.abs(Rcusp_b-Rc_m)
	idx = np.argmin(delta)
	prob_m = 1.0-float(idx)/Rcusp_b.size
	prob_sig = Rc_real/sig
	print prob,prob_p,prob_m
	print prob_sig


## fold
print "#### fold ####"

for i in range(len(p1)):
	pcen,Rf_real = p1[i],Rf[i]
	Rf_p,Rf_m = Rf[i]+Rfe[i],Rf[i]-Rfe[i]
	print pcen,Rf_real
	mask = np.logical_and(phi1>=pcen-2.5,phi1<pcen+2.5)

	Rfold_b = np.abs(Rfold[mask])
	sig = np.sqrt(Rfold_b.size)
	print Rfold_b.size, np.sqrt(Rfold_b.size)

	Rfold_b = np.sort(Rfold_b)
	delta = np.abs(Rfold_b-Rf_real)
	idx = np.argmin(delta)

	prob = 1.0-float(idx)/Rfold_b.size
	delta = np.abs(Rfold_b-Rf_p)
	idx = np.argmin(delta)
	prob_p = 1.0-float(idx)/Rfold_b.size

	delta = np.abs(Rfold_b-Rf_m)
	idx = np.argmin(delta)
	prob_m = 1.0-float(idx)/Rfold_b.size
	prob_sig = Rc_real/sig
	print prob,prob_p,prob_m
	print prob_sig


