import numpy as np
import commands
import matplotlib.pyplot as plt
import RcuspTool as rt

src_file = '/Volumes/sting_1/snap99_179899/179899_p1_src.dat'
src_rad = np.loadtxt(src_file)
src_px,src_py= np.degrees(src_rad[:,0])*3600,np.degrees(src_rad[:,1])*3600
src_g0 = np.array([ 5.946876e-02, -1.522927e-02])
src_p0 = np.array([src_px[0],src_py[0]])

src_tx,src_ty = (src_px-src_p0[0])+src_g0[0], (src_py-src_p0[1])+src_g0[1]

# SIE mag output
sie_out = open('/Volumes/sting_1/data/Rcusp_c/179899_p1_64_Rcusp_sie.txt','w')

for i in range(src_tx.size):
	findimg_file = open('179899_findimg.input','w')
	print src_tx[i],src_ty[i]
	findimg_file.write('gridmode 1\nset ngrid1=25\nset ngrid2=25\nset maxlev=3\n')
	findimg_file.write('startup 1 1\n')
	findimg_file.write('alpha 5.649761e-01 3.010972e-04 1.215803e-02 3.211276e-01 -7.083444e+01 1.016413e-01 -1.232670e+01 0.0 0.0 1.000000e+00\n')
	findimg_file.write('0 0 0 0 0 0 0 0 0 0\n')
	findimg_file.write('findimg '+str(src_tx[i])+' '+str(src_ty[i])+'\n')
	findimg_file.close()

	lens_cen = np.array([0.,0.])

	## get the output
		
	findimg_out = commands.getstatusoutput('./lensmodel 179899_findimg.input')
	lines = findimg_out[1].split('\n')
	last_line = lines[-2]
	quad_imghead = lines[-6]
	
	if(quad_imghead[0] == '#'):
		print quad_imghead
		raytrace_out = open('../../data/sub_gravlens/snap99_179899/sie_src'+str(i)+'.txt','w')
		raytrace_out.write('# '+str(src_tx[i])+'\t'+str(src_ty[i])+'\n')
		raytrace_out.write(quad_imghead+'\n')
		raytrace_out.write(lines[-5]+'\n'+lines[-4]+'\n'+lines[-3]+'\n'+lines[-2])
		raytrace_out.close()

		## read-in outputs and calculate Rcusp/Rfold
		find_table = np.loadtxt('../../data/sub_gravlens/snap99_179899/sie_src'+str(i)+'.txt',skiprows=1)

		img_x,img_y,img_f = find_table[:,0],find_table[:,1],find_table[:,2]
		
		rfold,rcusp,phi0,phi1 = rt.rcusp_tool(img_x,img_y,img_f,lens_cen)
		print rfold,rcusp,phi0,phi1
		sie_out.write(str(rfold)+'\t'+str(rcusp)+'\t'+str(phi0)+'\t'+str(phi1)+'\n')

		
		## find img B to assign parity

	else:
		print 'Not a quad'