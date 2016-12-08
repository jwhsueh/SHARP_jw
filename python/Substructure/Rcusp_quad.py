import pyfits
import numpy as np
import matplotlib.pyplot as plt
from photutils import find_peaks
from astropy.stats import sigma_clipped_stats
from itertools import combinations
from astropy import convolution

### change here ###
subID=263313
proj=1
NN='64'

n_src=100 # number of sources

filepath='/Volumes/sting_1/snap99_'+str(subID)
#output_name = '/Rcusp_'+str(subID)+'_p'+str(proj)+'_'+NN+'.txt' # no sub
output_name = '/Rcusp_'+str(subID)+'_p'+str(proj)+'sub_'+NN+'.txt' # w/ sub

rcusp_file=open(filepath+output_name,'w')
rcusp_file.write('# n_src = '+str(n_src)+'\n')
rcusp_file.write('# Rfold\tRcusp\n')

drop=0

for i in range(n_src):
	print "# "+str(i)+" source:"
	#filename='/image_'+str(subID)+'_p'+str(proj)+'_'+NN+'src_'+str(i)+'.fits' # no sub
	filename='/image_'+str(subID)+'_p'+str(proj)+'sub_'+NN+'src_'+str(i)+'.fits' # w/ sub

	image_fits=pyfits.open(filepath+filename)
	image=image_fits[0].data
	image_or=image

	## smoothing
	kernel=convolution.Gaussian2DKernel(3)
	image=convolution.convolve(image,kernel)

	# get image size
	img_size=image.shape[0]

	mean,median,std=sigma_clipped_stats(image,sigma=3.0)
	threshold=median+(20*std)
	peaks=find_peaks(image,threshold,box_size=5)

	print peaks

	#if (lens)
	# In the case that the fourth image is too small to see on the map
	# create an mock one at the furtherest corner
	if (len(peaks['x_peak']) == 3):
		print "The fourth image is too small to see on the map grid"
		corner_x = np.array([0,0,img_size,img_size])
		corner_y = np.array([0,img_size,0,img_size])

		# distance to four corners
		diff=np.zeros(4)
		for i in range(4):
			diff[i]=np.sum(np.abs(peaks['x_peak']-corner_x[i]))+np.sum(np.abs(peaks['y_peak']-corner_y[i]))

		idx = list(diff).index(max(diff))

		peaks.add_row([corner_x[idx],corner_y[idx],min(peaks['peak_value'])])


	lens_x,lens_y,lens_f=peaks['x_peak'],peaks['y_peak'],peaks['peak_value']

	if (len(lens_x)==4):

		## select image A,B & C [pick up the closet three pts]
		# find the shortest two line from C(4,2)
		quad_idx=np.arange(4)
		line_comb=list(combinations(quad_idx,2))
		line_len=np.zeros(len(line_comb))

		for i in range(len(line_comb)):
			pt0,pt1=line_comb[i][0],line_comb[i][1]
			line_len[i]=np.sqrt((lens_x[pt0]-lens_x[pt1])**2+(lens_y[pt0]-lens_y[pt1])**2)

		line_len_pop=list(line_len)
		min_0=np.min(line_len_pop) # 1st min, find the merging double
		idx0=list(line_len).index(min_0)

		line_len_pop.pop(idx0)
		min_1=np.min(line_len_pop) # 2nd min
		idx1=list(line_len).index(min_1)

		# get the line pts of 1st & 2nd mins

		line0,line1=line_comb[idx0],line_comb[idx1]
		tri_pt_idx=list(set().union(line0,line1))

		if (len(tri_pt_idx)>3):
			print "This is a CROSS lens!!"

		tri_mask=np.in1d(quad_idx,tri_pt_idx)
		tri_x,tri_y,tri_f=lens_x[tri_mask],lens_y[tri_mask],lens_f[tri_mask]
		## now we have the triple pts with smallest distance

		## save merging double index
		dou_mask=np.in1d(quad_idx,line0)
		dou_x,dou_y,dou_f=lens_x[dou_mask],lens_y[dou_mask],lens_f[dou_mask]

		## find image B (to assign parity)
		delta_x=np.max(tri_x)-np.min(tri_x)
		delta_y=np.max(tri_y)-np.min(tri_y)

		if (delta_x>delta_y): # projection axis=x
			proj_cord=tri_x
		else:
			proj_cord=tri_y

		middle=np.median(proj_cord)
		md_idx=list(proj_cord).index(middle)

		## assign parity
		mid_f=tri_f[md_idx] # this is the flux of img B (should assign negative parity)

		rfold=(np.sum(dou_f)-2.*mid_f)/np.sum(dou_f)
		rcusp=(np.sum(tri_f)-2.*mid_f)/np.sum(tri_f)

		print rfold,rcusp

		rcusp_file.write(str(rfold)+'\t'+str(rcusp)+'\n')

	else:
		print "Not detected as a quad!"
		drop=drop+1
		rcusp_file.write('# '+str(i)+' source dropped\n')

	#if(len(peaks['x_peak'])!=3 & len(peaks['x_peak'])!=4):
	#	print "Not a quad"


print '# of drop ='+str(drop)

