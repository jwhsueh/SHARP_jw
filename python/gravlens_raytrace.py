import numpy as np
import commands

model_path = '/Users/jwhsueh/Documents/SHARP_jw/models/snap99_285514/'

## construct source array
xx = np.linspace(110,160,51)
yy = np.linspace(110,150,41)
sx,sy = np.meshgrid(xx,yy)

sx,sy=sx.flatten(),sy.flatten()

findimg_file = open('285514_findimg.input','w')

for i in range(sx.size):

	findimg_file.write('start up 1 1\n')
	findimg_file.write('alpha 5.463769e+01 1.366652e+02 1.259495e+02 4.009148e-01 6.495593e+01 0.0 0.0 0.0 0.0 1.000000e+00 \n')
	findimg_file.write('0 0 0 0 0 0 0 0 0 0\n')
	findimg_file.write('findimg '+str(sx[i])+' '+str(sy[i])+'\n')

	## get the output

	findimg_out = commands.getstatusoutput('lensmodel 285514_findimg.input')
	print findimg_out
