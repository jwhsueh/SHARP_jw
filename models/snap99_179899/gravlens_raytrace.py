import numpy as np
import commands
import matplotlib.pyplot as plt
import RcuspTool as rt

model_path = '/Users/jwhsueh/Documents/SHARP_jw/models/snap99_179899/'

## construct source array
res = 51
xx = np.linspace(0.7,1.0,res)
yy = np.linspace(0.7,1.0,res)
#sx,sy = np.meshgrid(xx,yy)
mag_3 = np.zeros((res,res))
lens_cen = np.array([8.331855e-01, 8.614601e-01])

#sx,sy=sx.flatten(),sy.flatten()
rc,rf,p0,p1 = np.zeros((res,res)),np.zeros((res,res)),np.zeros((res,res)),np.zeros((res,res))

src_file = '/Volumes/sting_1/snap99_179899/179899_p1_src.dat'
src_rad = np.loadtxt(src_file)
src_px,src_py= np.degrees(src_rad[:,0])*3600,np.degrees(src_rad[:,1])*3600
src_g0 = np.array([ 5.946876e-02, -1.522927e-02])
src_p0 = np.array([src_px[0],src_py[0]])

src_tx,src_ty = (src_px-src_p0[0])+src_g0[0], (src_py-src_p0[1])+src_g0[1]

src_out = open('src_pt.txt','w')
for i in range(src_ty.size):
	src_out.write(str(src_tx[i])+'\t'+str(src_ty[i])+'\n')

n_src = 0

for i in range(xx.size):
	for j in range(yy.size):

		findimg_file = open('179899_findimg.input','w')
		print xx[i],yy[j]
		findimg_file.write('gridmode 1\nset ngrid1=25\nset ngrid2=25\nset maxlev=3\n')
		findimg_file.write('startup 1 1\n')
		findimg_file.write('alpha 5.041207e-01 8.331855e-01 8.614601e-01 3.763852e-01 -6.631899e+01 9.905569e-02 -8.770241e+00 0.0 0.0 1.000000e+00 \n')
		findimg_file.write('0 0 0 0 0 0 0 0 0 0\n')
		findimg_file.write('findimg '+str(xx[i])+' '+str(yy[j])+'\n')
		findimg_file.close()

		## get the output
		
		findimg_out = commands.getstatusoutput('./lensmodel 179899_findimg.input')
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

print n_src

np.savetxt('179899_rfold.txt',rf)
np.savetxt('179899_rcusp.txt',rc)
np.savetxt('179899_phi0.txt',p0)
np.savetxt('179899_phi1.txt',p1)


#plt.
plt.imshow(rf,extent=[xx[0],xx[-1],yy[0],yy[-1]])
plt.scatter(src_tx,src_ty,marker='x',s=1)
plt.colorbar()
plt.show()