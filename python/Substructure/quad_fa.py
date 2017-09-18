from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from photutils import find_peaks
from astropy.stats import sigma_clipped_stats
from itertools import combinations
from astropy.modeling import models, fitting 
from astropy import convolution
from astropy.table import Table

n_real = 1
n_start = 0

re = 0.75
N_pix = 512
pix_size = 0.75*2*2/N_pix #arcsec

mer_beam = 50 # mas
beam = mer_beam/(pix_size*1000)

print beam

#imgpath='/Volumes/sting_1/subs/glamer_out/B1422_0010/'
imgpath='/Volumes/sting_1/subs/'
outfile = imgpath+'B1422_0010_glamer_chain.txt'
#output = open(outfile,'w')

img_idx = 2
f_idx = 3
boxsize = 20

for i in range(n_real):
	realID = i+n_start

	#image_fits=fits.open(imgpath+'image_fsub1_'+str(realID)+'.fits')
	image_fits=fits.open(imgpath+'image_sie_'+str(realID)+'.fits')
	image=image_fits[0].data

	## smoothing
	kernel=convolution.Gaussian2DKernel(beam)
	image=convolution.convolve(image,kernel)

	mean,median,std=sigma_clipped_stats(image,sigma=3.0)
	threshold=median+(20*std)
	peaks=find_peaks(image,threshold,box_size=5) # w/ convolve 

	#print peaks

	x_pix,y_pix = np.array(peaks['x_peak']),np.array(peaks['y_peak'])
	f_img = np.array(peaks['peak_value'])

	if (len(x_pix)==4):
		x_img,y_img = x_pix-x_pix[img_idx],y_pix-y_pix[img_idx]
		

		## ------ Gaussian fit	
		## img cut
		gaus2d = models.Gaussian2D(amplitude=1., x_mean=0, y_mean=0., x_stddev=1.,y_stddev=1.)
		fit_g = fitting.LevMarLSQFitter()
		
		grid_y,grid_x = np.mgrid[:boxsize,:boxsize]
		#grid_y,grid_x = np.mgrid[:img_size,:img_size]
		grid_yy,grid_xx = np.mgrid[-2*boxsize:2*boxsize,-2*boxsize:2*boxsize]

		for k in range(4):
			img_cut = image[y_pix[k]-boxsize/2:y_pix[k]+boxsize/2,x_pix[k]-boxsize/2:x_pix[k]+boxsize/2]
			fit_out = fit_g(gaus2d,grid_x,grid_y,img_cut)
			#max_v = np.max(fit_out(grid_x,grid_y))
			fit_mod = fit_out(grid_x,grid_y)
						
			#plt.imshow(fit_out(grid_xx,grid_yy))
			#plt.colorbar()
			#plt.show()
			#plt.clf()
			
			ga_array = fit_out(grid_xx,grid_yy)
			hm = ga_array[ga_array>np.max(ga_array)/2]
			f_img[k] = np.sum(hm)
		## gauss fit -- end

		f_img = f_img/f_img[f_idx]
		x_img,y_img = x_img*pix_size,y_img*pix_size




		print x_img,y_img
		print f_img

		#output.write(str(x_img[3])+'\t'+str(y_img[3])+'\t'+str(f_img[2])+'\t'+str(x_img[1])+'\t'+str(y_img[1])+'\t'+str(f_img[1])+'\t'+str(x_img[0])+'\t'+str(y_img[0])+'\t'+str(f_img[0])+'\n')

	else:
		print realID,len(x_img)


