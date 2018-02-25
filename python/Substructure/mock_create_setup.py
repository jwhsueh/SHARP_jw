import numpy as np
import gravlens_tool as gt
import commands

N=90 # number of mocks
nint=10
path = '/Volumes/sting_1/subs/mock_test/'

## macro model variable setup
## gaussaian
q = np.random.rand(N)*0.2+0.25 # elliplicity flat dist [0.25,0.45]
gamma = np.random.normal(0.05,0.01,N)
th_gamma = np.random.rand(N)*180
b = np.random.rand(N)*0.4+0.6 # Einstein radius flat distribution [0.6,1.0]

## setup gravlens macro model array
m1_table = np.zeros((N,10))
m1_table[:,0] = b
m1_table[:,3] = q
m1_table[:,5] = gamma
m1_table[:,6] = th_gamma
m1_table[:,9] = 1.0

np.savetxt(path+'m1_table2.txt',m1_table,fmt='%f')

## mock lens redshift

zl = np.random.rand(N)*0.3+0.3
zs = np.random.rand(N)+1.0

np.savetxt(path+'m1_redshift2.txt',np.c_[zl,zs],fmt='%f')

## ----- create gravlens files to create caustics ----- ##
path1 = path+'gravlens_file/'
for i in range(N):
	macro_mod = m1_table[i,:]
	z = np.array([zl[i],zs[i]])
	midx = i+nint
	gt.create_caustic(macro_mod,z,'plotcrit.input','mock'+str(midx)+'_crit.dat',path1)
	opt_out = commands.getstatusoutput('./lensmodel '+path1+'plotcrit.input')
	#rint opt_out

## ----- create src file ----- ##

src_array = np.zeros((N,2))
for i in range(N):
	midx = i+nint
	sx,sy = gt.quad_src(1,np.array([0.0,0.0]),path1,'mock'+str(midx)+'_crit.dat') #nsrc, lens_center

	src_array[i,0] = sx[0]
	src_array[i,1] = sy[0]

np.savetxt(path+'m1_src2.txt',src_array,fmt='%f')

	