import numpy as np
import commands
import matplotlib.pyplot as plt
import RcuspTool as rt

model_path = '/Users/jwhsueh/Documents/SHARP_jw/models/snap99_179899/'
lens_path = '/Volumes/sting_1/snap99_179899/'

sublist='179899_p1'
subid='179899'
n_src = 100

rcusp_file = open(model_path+sublist+'_sie_rcusp_sh.txt','w')
ast_file = open(model_path+sublist+'_sie_ast_sh.txt','w')

for i in range(n_src):
	print '#### nrc '+str(i)+' ####'

	## opt src file
	imgfile = lens_path+'image_'+sublist+'_64src_'+str(i)+'_s_wcs.txt'
	tab=np.genfromtxt(imgfile)
	print tab
	opt_src = open(model_path+subid+'_src.dat','w')

	opt_src.write('1\n')
	opt_src.write('0 0 100000\n')
	opt_src.write('0.3 100000\n')
	opt_src.write('99 100000\n')
	opt_src.write('0.2 100000\n\n')
	opt_src.write('1\n4\n')
	opt_src.write(str(tab[0])[1:-1]+' 1.368  0.006 1000  0.0 0.0 \n')
	opt_src.write(str(tab[1])[1:-1]+' 1.368  0.006 1000  0.0 0.0 \n')
	opt_src.write(str(tab[2])[1:-1]+' 1.368  0.006 1000  0.0 0.0 \n')
	opt_src.write(str(tab[3])[1:-1]+' 1.368  0.006 1000  0.0 0.0 \n')

	opt_src.close()

	## run opt 
	
	opt_out = commands.getstatusoutput('./lensmodel opt_179899_gravlens.input')
	opt_out = opt_out[1].split('\n')
	print opt_out[-5]

	## read best.dat
	# img
	bestfile = open(model_path+'best.dat','r')
	best_line = []
	for line in bestfile.readlines():
		best_line.append(line)
	##

	lens_model = best_line[1].split()
	print lens_model
	shear = float(lens_model[6])
	len_x,len_y = float(lens_model[2]),float(lens_model[3])
	len_cen = np.array([len_x,len_y])
	print 'shear= '+str(shear)

	# chi2

	chi2 = best_line[4].split()
	chi2 =float(chi2[2])
	print 'chi-square='+str(chi2)

	img_x,img_y,img_f=np.empty(4),np.empty(4),np.empty(4)
	ray_x,ray_y=np.empty(4),np.empty(4)

	for j in range(4):
		find_line = best_line[-3-j].split(' ')
		img_line = []
		for element in find_line:
			if element!='':
				img_line.append(element)

		#print img_line
		img_x[j] = float(img_line[-4])
		img_y[j] = float(img_line[-3])
		img_f[j] = float(img_line[-2])

		ray_x[j] = float(img_line[0])
		ray_y[j] = float(img_line[1])

	print img_x

	## astrometry anomaly
	ast = np.sum(np.sqrt((img_x-ray_x)**2+(img_y-ray_y)**2))
	ast_file.write(str(ast)+'\t'+str(shear)+'\t'+str(chi2)+'\n')
	print ast
	
	## Rcusp calculation
	
	rfold,rcusp,phi0,phi1 = rt.rcusp_tool(img_x,img_y,img_f,len_cen)
	rcusp_file.write(str(rfold)+'\t'+str(rcusp)+'\t'+str(phi0)+'\t'+str(phi1)+'\n')
	print rfold,rcusp,phi0,phi1

	n_src = n_src+1

ast_file.close()
rcusp_file.close()
