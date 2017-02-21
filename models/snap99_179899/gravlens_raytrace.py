import numpy as np
import commands
import matplotlib.pyplot as plt

model_path = '/Users/jwhsueh/Documents/SHARP_jw/models/snap99_179899/'

## construct source array
res = 51
xx = np.linspace(0.7,1.0,res)
yy = np.linspace(0.7,1.0,res)
#sx,sy = np.meshgrid(xx,yy)
mag_3 = np.zeros((res,res))

#sx,sy=sx.flatten(),sy.flatten()

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
			raytrace_out = open('./rt_gravlens/rt_src'+str(n_src)+'.txt','w')
			raytrace_out.write(str(xx[i])+'\t'+str(yy[j])+'\n')
			raytrace_out.write(quad_imghead+'\n')
			raytrace_out.write(lines[-5]+'\n'+lines[-4]+'\n'+lines[-3]+'\n'+lines[-2])
			raytrace_out.close()

			## read-in outputs and generate magnification on grid
			table = np.loadtxt('./rt_gravlens/rt_src'+str(n_src)+'.txt',skiprows=1)

			mag = table[:,2]
			mag = list(mag)
			mag.pop(mag.index(np.min(mag)))
			mag_3[j,i] = np.min(mag)

			n_src = n_src+1

		else:
			print 'Not a quad'

print n_src

np.savetxt('./rt_gravlens/mag_3.txt',mag_3)


#plt.
plt.imshow(mag_3)
plt.colorbar()
plt.show()