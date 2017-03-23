import numpy as np
import matplotlib.pyplot as plt
import commands
import RcuspTool as rt

tab = np.loadtxt('1798992_crit.dat')
ux,uy = np.append(tab[:,2],tab[:,6]),np.append(tab[:,3],tab[:,7])

lens_cen = np.array([8.331855e-01, 8.614601e-01])

ur = np.sqrt((ux-lens_cen[0])**2+(uy-lens_cen[1])**2)
mask = ur<0.2
ur = ur[mask]

rcusp_out = open('1798992_area15_rcusp.txt','w')
src_out = open('1798992_area15_src.txt','w')

ux,uy = ux[mask],uy[mask]

#plt.scatter(ux,uy,s=1)
#plt.show()

## calculate r_caus(min)
r_min = np.min(ur)

## draw randomly
n_src = 500
i=0
sx,sy = [],[]
sx0,sy0 = 0.0,0.0
while (i<n_src):
	xi = np.random.random_sample()*(0.3)+0.7
	yi = np.random.random_sample()*(0.3)+0.7

	# distance to caustic
	ri = np.sqrt((xi-ux)**2+(yi-uy)**2)
	rci = np.sqrt((xi-lens_cen[0])**2+(yi-lens_cen[1])**2)
	dist = np.min(ri)
	idx = np.argmin(ri)
	rc = ur[idx]

	if np.logical_and(dist<0.15*r_min,dist>0.01*r_min):
		if (rci<rc):
			if np.logical_and(xi!=sx0,yi!=sy0):
				sx.append(xi)
				sy.append(yi)

				sx0,sy0=xi,yi
				i=i+1

sx,sy = np.array(sx),np.array(sy)
src_out.write('# '+str(n_src)+' src\n')
for i in range(sx.size):
	src_out.write(str(sx)+'\t'+str(sy)+'\n')

## gravlens findimg


for i in range(sx.size):

	findimg_file = open('1798992_findimg.input','w')
	findimg_file.write('gridmode 1\nset ngrid1=25\nset ngrid2=25\nset maxlev=3\n')
	findimg_file.write('startup 1 1\n')
	findimg_file.write('alpha 5.041207e-01 8.331855e-01 8.614601e-01 1.763852e-01 -6.631899e+01 0.0 0.0 0.0 0.0 1.000000e+00 \n')
	#findimg_file.write('expdiskA 3.341240e-01 8.331855e-01 8.614601e-01 6.581240e-01 -5.995461e+01 0.0 0.0 4.114300e-01 0.0 0.0 \n')
	findimg_file.write('0 0 0 0 0 0 0 0 0 0\n')
	#findimg_file.write('0 0 0 0 0 0 0 0 0 0\n')
	findimg_file.write('findimg '+str(sx[i])+' '+str(sy[i])+'\n')
	findimg_file.close()

	## get the output
	
	findimg_out = commands.getstatusoutput('./lensmodel 1798992_findimg.input')
	lines = findimg_out[1].split('\n')
	last_line = lines[-2]
	quad_imghead = lines[-6]
	
	if(quad_imghead[0] == '#'):
		print quad_imghead
		raytrace_out = open('../../data/sub_gravlens/snap99_1798992/area_15/a15_src'+str(i)+'.txt','w')
		raytrace_out.write(str(sx[i])+'\t'+str(sy[i])+'\n')
		raytrace_out.write(quad_imghead+'\n')
		raytrace_out.write(lines[-5]+'\n'+lines[-4]+'\n'+lines[-3]+'\n'+lines[-2])
		raytrace_out.close()

		## read-in outputs and generate magnification on grid
		tab = np.loadtxt('../../data/sub_gravlens/snap99_1798992/area_15/a15_src'+str(i)+'.txt',skiprows=1)
		print tab
		img_x,img_y,img_f = tab[:,0],tab[:,1],tab[:,2]

		rfold,rcusp,phi0,phi1 = rt.rcusp_tool(img_x,img_y,img_f,lens_cen)
		#rc[i,j],rf[i,j],p0[i,j],p1[i,j] = rcusp,rfold,phi0,phi1

		rcusp_out.write(str(rfold)+'\t'+str(rcusp)+'\t'+str(phi0)+'\t'+str(phi1)+'\n')

	else:
		rcusp_out.write('# '+str(i)+' drop \n')
		print 'Not a quad'

rcusp_out.close()
'''
tab = np.loadtxt('1798992_area15_rcusp.txt')
plt.scatter(tab[:,2],np.abs(tab[:,1]),color='b',marker='*',label='e = 0.18')
plt.ylim(0,0.7)
plt.xlim(0,200)
plt.title('edge-on Rcusp 0.35*r')
plt.ylabel('delta phi')
plt.xlabel('|Rcusp|')
plt.legend(loc=2)
#plt.show()
plt.savefig('179899_d3_a15_e66_rcusp.png')
'''
