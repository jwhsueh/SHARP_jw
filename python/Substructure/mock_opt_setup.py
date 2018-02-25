import numpy as np
import gravlens_tool as gt


nint=0
N=100 - nint # number of mocks
path = '/Volumes/sting_1/subs/mock_test/'
path1 = path+'gravlens_file/'

m1_model =  np.loadtxt(path+'m1_table.txt')
m1_obs = np.genfromtxt(path+'m1_raytrace.txt',dtype=float)

## m2 files
m2_model = open(path+'m2_table.txt','w')
m2_src = open(path+'m2_src.txt','w')

drop = np.loadtxt(path+'m2_dropout.txt')
drop_model = np.zeros(7)
drop_src = np.zeros(2)

for i in range(N):
	midx = i+nint
	print '#######'
	print 'mock '+str(midx)
	print '#######'

	if np.in1d(midx,drop):
		m2_model.write(str(drop_model[0])+'\t'+str(drop_model[1])+'\t'+str(drop_model[2])+'\t'+str(drop_model[3])+'\t'+str(drop_model[4])+'\t'+str(drop_model[5])+'\t'+str(drop_model[6])+'\n')
		m2_src.write(str(drop_src)[1:-1]+'\n')

	else:
		## drop out file
		mock_drop = open(path+'m2_dropout.txt','a')

		macro_model = m1_model[midx,:]
		mock_obs = m1_obs[midx,:]

		img_x,img_y,img_f = mock_obs[:4],mock_obs[4:8],mock_obs[8:12]
		img_f = np.abs(img_f)
		err_x = 0.0003
		#err_f = img_f*0.03

		## create mock dat
		dat_file = 'mock'+str(midx)+'_obs.dat'
		gt.create_mockdat(img_x,img_y,img_f,err_x,path1,dat_file)

		## create opt file
		opt_file = 'opt_mocks.input'
		gt.create_sieopt(macro_model,path1,dat_file,opt_file)

		## do optimization
		best_line,lens_para,src,chi2 = gt.gravlens_opt(path1,opt_file)
		if chi2>1.0:
			mock_drop.write(str(midx)+'\n')

		m2_model.write(str(lens_para[0])+'\t'+str(lens_para[1])+'\t'+str(lens_para[2])+'\t'+str(lens_para[3])+'\t'+str(lens_para[4])+'\t'+str(lens_para[5])+'\t'+str(lens_para[6])+'\n')
		m2_src.write(str(src)[1:-1]+'\n')

		mock_drop.close()

m2_model.close()
m2_src.close()
