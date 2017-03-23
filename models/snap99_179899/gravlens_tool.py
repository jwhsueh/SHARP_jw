import numpy as np
import commands

def findimg_quad(lens,src,target):
	# lens is the line describe lens setup in gravlens
	# src is an np array record the source position
	# target is the name of the lens

	filename = target+'_findimg.input'
	findimg_file = open(filename,'w')

	findimg_file.write('gridmode 1\nset ngrid1=25\nset ngrid2=25\nset maxlev=3\n')
	findimg_file.write('startup 1 1\n')
	findimg_file.write(lens+'\n')
	findimg_file.write('0 0 0 0 0 0 0 0 0 0\n')
	findimg_file.write('findimg '+str(src[0])+' '+str(src[1])+'\n')
	findimg_file.close()

	## get output
	findimg_out = commands.getstatusoutput('./lensmodel '+filename)
	lines = findimg_out[1].split('\n')
	last_line = lines[-2]
	quad_imghead = lines[-6]

	if(quad_imghead[0] == '#'):
			print quad_imghead
			raytrace_out = open('../../data/sub_gravlens/snap99_179899/rt_src'+str(n_src)+'.txt','w')
			raytrace_out.write(str(xx[i])+'\t'+str(yy[j])+'\n')
			raytrace_out.write(quad_imghead+'\n')
			raytrace_out.write(lines[-5]+'\n'+lines[-4]+'\n'+lines[-3]+'\n'+lines[-2])
			raytrace_out.close()

			## read-in outputs and generate magnification on grid
			tab = np.loadtxt('../../data/sub_gravlens/snap99_179899/rt_src'+str(n_src)+'.txt',skiprows=1)
			print tab
			img_x,img_y,img_f = tab[:,0],tab[:,1],tab[:,2]

			rfold,rcusp,phi0,phi1 = rt.rcusp_tool(img_x,img_y,img_f,lens_cen)
			rc[i,j],rf[i,j],p0[i,j],p1[i,j] = rcusp,rfold,phi0,phi1

			n_src = n_src+1

		else:
			print 'Not a quad'