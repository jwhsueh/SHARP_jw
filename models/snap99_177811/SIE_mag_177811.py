import numpy as np
import commands
import matplotlib.pyplot as plt
import RcuspTool as rt

src_file = '/Volumes/sting_1/snap99_177811/177811_p2_src.dat'
src_rad = np.loadtxt(src_file)
src_px,src_py= np.degrees(src_rad[:,0])*3600,np.degrees(src_rad[:,1])*3600
src_g0 = np.array([ -1.862299e-02, -1.636522e-02])
src_p0 = np.array([src_px[0],src_py[0]])

src_tx,src_ty = (src_px-src_p0[0])+src_g0[0], (src_py-src_p0[1])+src_g0[1]

# SIE mag output
sie_out = open('/Volumes/sting_1/data/Rcusp_c/177811_p2_64_Rcusp_sie.txt','w')

for i in range(src_tx.size):
	findimg_file = open('177811_findimg.input','w')
	print src_tx[i],src_ty[i]
	findimg_file.write('gridmode 1\nset ngrid1=25\nset ngrid2=25\nset maxlev=3\n')
	findimg_file.write('startup 1 1\n')
	findimg_file.write('alpha 4.110707e-01 6.231697e-02 3.422431e-03 5.558948e-01 5.618475e+01 8.076018e-02 -1.307103e+01 0.0 0.0 1.000000e+00\n')
	findimg_file.write('0 0 0 0 0 0 0 0 0 0\n')
	findimg_file.write('findimg '+str(src_tx[i])+' '+str(src_ty[i])+'\n')
	findimg_file.close()

	lens_cen = np.array([0.,0.])

	## get the output
		
	findimg_out = commands.getstatusoutput('./lensmodel 177811_findimg.input')
	lines = findimg_out[1].split('\n')
	#print lines
	last_line = lines[-2]
	quad_imghead = lines[-6]
	
	if(quad_imghead[0] == '#'):
		print quad_imghead
		raytrace_out = open('../../data/sub_gravlens/snap99_177811/sie_src'+str(i)+'.txt','w')
		raytrace_out.write('# '+str(src_tx[i])+'\t'+str(src_ty[i])+'\n')
		raytrace_out.write(quad_imghead+'\n')
		raytrace_out.write(lines[-5]+'\n'+lines[-4]+'\n'+lines[-3]+'\n'+lines[-2])
		raytrace_out.close()

		## read-in outputs and calculate Rcusp/Rfold
		find_table = np.loadtxt('../../data/sub_gravlens/snap99_177811/sie_src'+str(i)+'.txt',skiprows=1)

		img_x,img_y,img_f = find_table[:,0],find_table[:,1],find_table[:,2]
		
		rfold,rcusp,phi0,phi1 = rt.rcusp_tool(img_x,img_y,img_f,lens_cen)
		print rfold,rcusp,phi0,phi1
		sie_out.write(str(rfold)+'\t'+str(rcusp)+'\t'+str(phi0)+'\t'+str(phi1)+'\n')

		
		## find img B to assign parity

	else:
		print 'Not a quad'