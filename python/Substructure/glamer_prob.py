import numpy as np

path='/Volumes/sting_1/data/'

#list_file = path+'snap99_elp2_Rcusp.txt'
list_file = path+'snap99_sn2_Rcusp.txt'
sie_path='/Volumes/sting_1/snap99_result/snap99_111/'
sie_path2='/Volumes/sting_1/snap99_result/snap99_222/'

file_list = np.genfromtxt(list_file,dtype='str')
#file_list2 = np.genfromtxt(list_file2,dtype='str')

Rfold,Rcusp,phi0,phi1 = np.empty(1),np.empty(1),np.empty(1),np.empty(1)
for file_one in file_list:
	#table = np.loadtxt(path+'Rcusp_f/'+file_one+'_64_Rcusp_ga.txt')
	table = np.loadtxt(sie_path+file_one+'_Rcusp_ga.txt') # sie
	Rfold,Rcusp,phi0,phi1 = np.append(Rfold,table[:,0]),np.append(Rcusp,table[:,1]),np.append(phi0,table[:,2]),np.append(phi1,table[:,3])
'''
for file_one in file_list2:
	#table = np.loadtxt(path+'Rcusp_f/'+file_one+'_64_Rcusp_ga.txt')
	table = np.loadtxt(sie_path2+file_one+'_Rcusp_ga.txt') # sie
	Rfold,Rcusp,phi0,phi1 = np.append(Rfold,table[:,0]),np.append(Rcusp,table[:,1]),np.append(phi0,table[:,2]),np.append(phi1,table[:,3])
'''

mask = np.abs(Rfold)<0.7
Rfold,Rcusp,phi0,phi1=Rfold[mask],Rcusp[mask],phi0[mask],phi1[mask]
mask = np.abs(Rcusp)<0.7
Rfold,Rcusp,phi0,phi1=Rfold[mask],Rcusp[mask],phi0[mask],phi1[mask]

#Rfold[phi1<27] = Rfold[phi1<27]/2

## ---- real lens data
# elp
tab_e = np.loadtxt('../../data/flux_ratio.txt')
Rf_e,Rc_e,p0_e,p1_e,Rfe_e,Rce_e = tab_e[:,3],tab_e[:,1],tab_e[:,0],tab_e[:,2],tab_e[:,4],tab_e[:,5]

tab_d = np.loadtxt('../../data/flux_ratio_disk.txt')
Rf_d,Rc_d,p0_d,p1_d,Rfe_d,Rce_d = tab_d[:,3],tab_d[:,1],tab_d[:,0],tab_d[:,2],tab_d[:,4],tab_d[:,5]

## change here

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
	mask_m,mask_p = Rcusp_b>=Rc_m,Rcusp_b>=Rc_p
	Rcusp_m,Rcusp_p = Rcusp_b[mask_m],Rcusp_b[mask_p]
	print Rcusp_b.size
	sig = np.sqrt(Rcusp_b.size)
	delta_sig = np.sqrt(Rcusp_m.size-Rcusp_p.size)
	#sig_m,sig_p = np.sqrt(Rcusp_m.size),np.sqrt(Rcusp_p.size)

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
	#poi_m,poi_p = Rc_real/sig_m,Rc_real/sig_p
	print prob,prob_p,prob_m
	print prob_sig
	print "poison noise"
	print Rc_real/delta_sig


## fold
print "#### fold ####"

for i in range(len(p1)):
	pcen,Rf_real = p1[i],Rf[i]
	Rf_p,Rf_m = Rf[i]+Rfe[i],Rf[i]-Rfe[i]
	print pcen,Rf_real
	mask = np.logical_and(phi1>=pcen-2.5,phi1<pcen+2.5)

	Rfold_b = np.abs(Rfold[mask])
	sig = np.sqrt(Rfold_b.size)
	#print Rfold_b.size, np.sqrt(Rfold_b.size)
	mask_m,mask_p = Rfold_b>=Rf_m,Rfold_b>=Rf_p
	Rfold_m,Rfold_p = Rfold_b[mask_m],Rfold_b[mask_p]
	sig_m,sig_p = np.sqrt(Rfold_m.size),np.sqrt(Rfold_p.size)
	delta_sig = np.sqrt(Rfold_m.size-Rfold_p.size)

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
	prob_sig = Rf_real/sig
	poi_m,poi_p = Rf_real/sig_m,Rf_real/sig_p
	print prob,prob_p,prob_m
	print prob_sig
	print "poison noise"
	print Rf_real/delta_sig


