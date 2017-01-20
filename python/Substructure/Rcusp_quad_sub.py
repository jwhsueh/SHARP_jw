import pyfits
import numpy as np
import matplotlib.pyplot as plt
from photutils import find_peaks
from astropy.stats import sigma_clipped_stats
from itertools import combinations
from astropy import convolution
from astropy.table import Table

### change here ###
subID=238069
proj=2
NN='64'

rmax = 4.14821e-06
real_size = 1.1*2*rmax # rad
real_size = np.degrees(real_size)*3600*1000 # mas
print real_size
mer_beam = 10 # mas

n_src=50 # number of sources
img_size=256

beam = mer_beam/(real_size/img_size)
print beam

filepath='/Volumes/sting_1/snap99_'+str(subID)
output_name = '/'+str(subID)+'_p'+str(proj)+'_'+NN+'_Rcusp.txt' # no sub
#output_name = '/'+str(subID)+'_p'+str(proj)+'sub_'+NN+'_Rcusp.txt' # w/ sub

rcusp_file=open(filepath+output_name,'w')
rcusp_file.write('# n_src = '+str(n_src)+'\n')
rcusp_file.write('# Rfold\tRcusp\tRfold_s\tRcusp_s\n')


drop=0
drop_file='/'+str(subID)+'_p'+str(proj)+'_'+NN+'_drop.txt'
#drop_log=open(filepath+drop_file,'w')
drop_log=open(filepath+drop_file,'a+')

## suscess flag
nosub=0
wsub=0  # 1 = drop


for i in range(n_src):
	print "# "+str(i)+" source:"

	# go over nosub
	filename='/image_'+str(subID)+'_p'+str(proj)+'_'+NN+'src_'+str(i)+'.fits' # no sub
	filename_sub='/image_'+str(subID)+'_p'+str(proj)+'sub_'+NN+'src_'+str(i)+'.fits' # w/ sub

	image_fits=pyfits.open(filepath+filename)
	image=image_fits[0].data
	image_or=image

	## smoothing
	kernel=convolution.Gaussian2DKernel(beam)
	image=convolution.convolve(image,kernel)

	mean,median,std=sigma_clipped_stats(image,sigma=3.0)
	threshold=median+(20*std)
	peaks=find_peaks(image,threshold,box_size=5)

	lens_x,lens_y,lens_f = peaks['x_peak'],peaks['y_peak'],peaks['peak_value']
	#print peaks

	# In case that the fourth image is too small to see on the map
	# create an mock one at the furtherest corner
	if (len(peaks['x_peak']) == 3):
		print "The fourth image is too small to see on the map grid"
		corner_x = np.array([0,0,img_size-1,img_size-1])
		corner_y = np.array([0,img_size-1,0,img_size-1])

		# distance to four corners
		diff=np.zeros(4)
		for j in range(4):
			diff[j]=np.sum(np.abs(peaks['x_peak']-corner_x[j]))+np.sum(np.abs(peaks['y_peak']-corner_y[j]))

		idx = list(diff).index(max(diff))

		peaks.add_row([corner_x[idx],corner_y[idx],min(peaks['peak_value'])])
		lens_x,lens_y,lens_f=peaks['x_peak'],peaks['y_peak'],peaks['peak_value']


	# In case that source sometimes show up on the grid
	if (len(peaks['x_peak']) == 5):
		print "mask out the source"
		mask = list(peaks['peak_value']).index(min(peaks['peak_value']))
		peaks.remove_row(mask)
		lens_x,lens_y,lens_f=peaks['x_peak'],peaks['y_peak'],peaks['peak_value']
	
	if (len(peaks['x_peak']) >4):
		print "Go into clustering"
		## ----- Clustering process ----
		#	This part pick C(n,4) pts from peak detection result
		#	active pts: 4 pts that has been selected
		# 	rest pts: pts other than 4 active pts
		#	calculate the minimum distance to any active pts that is "brighter than" this rest pt
		#	If there's no active pts brighter than this rest pt, abandon this combination
		#	Survival combination is the one has minimum sum of clustering distance (w/ flux weighted)
		#	flux wighted: devided by the total flux of active pts. this means when two sets have close minimum sum of clustering
		#				the set with larger total flux will survive.
		## ----------------------

		img_idx=np.arange(len(peaks['x_peak']))
		comb_list=list(combinations(img_idx,4))
		#print comb_list
		comb_dist_sum=np.zeros(len(comb_list))	# 1e10: this combination is abandoned 

		for k in range(len(comb_list)):
			active_idx=comb_list[k]
			act_mask=np.in1d(img_idx,active_idx)
			rest_idx=img_idx[~act_mask]

			#print active_idx,rest_idx
			act_x,act_y,act_f=peaks['x_peak'][act_mask],peaks['y_peak'][act_mask],peaks['peak_value'][act_mask]
			rest_x,rest_y,rest_f=peaks['x_peak'][~act_mask],peaks['y_peak'][~act_mask],peaks['peak_value'][~act_mask]

			## start to find brighter points than each rest points
			# compare maximum peak in active and rest pts pool
			
			if(np.max(rest_f)>np.max(act_f)): # rest pts actually contain one of lens image [This is not enough!]
				comb_dist_sum[k]=1e10
				# this set of active pts is abandoned

			else:
				## ---- start to calculate minimum dist to one of the active pts brighter than this rest pt

				dist_sum=0.

				for j in range(len(rest_x)):
					# who is brighter
					flux_mask=act_f>rest_f[j]
					br_x,br_y=act_x[flux_mask],act_y[flux_mask]

					dist_to_bpt=np.sqrt((rest_x[j]-br_x)**2+(rest_y[j]-br_y)**2)
					dist_sum=dist_sum+np.min(dist_to_bpt)

				comb_dist_sum[k]=dist_sum/np.sum(act_f)
			
			#print np.max(rest_f),np.max(act_f)
			#print comb_dist_sum[i]

		sur_comb_idx = list(comb_dist_sum).index(np.min(comb_dist_sum))
		sur_mask = np.in1d(img_idx,comb_list[sur_comb_idx])

		lens_x,lens_y,lens_f=peaks['x_peak'][sur_mask],peaks['y_peak'][sur_mask],peaks['peak_value'][sur_mask]

	### end of clustering ###

	quad_table=Table([lens_x,lens_y,lens_f])
	print quad_table
	

	if (len(lens_x)==4):

		## select image A,B & C [pick up the closet three pts]
		# find the shortest two line from C(4,2)
		quad_idx=np.arange(4)
		line_comb=list(combinations(quad_idx,2))
		line_len=np.zeros(len(line_comb))

		for j in range(len(line_comb)):
			pt0,pt1=line_comb[j][0],line_comb[j][1]
			line_len[j]=np.sqrt((lens_x[pt0]-lens_x[pt1])**2+(lens_y[pt0]-lens_y[pt1])**2)

		line_len_pop=list(line_len)
		min_0=np.min(line_len_pop) # 1st min, find the merging double
		idx0=list(line_len).index(min_0)

		line_len_pop.pop(idx0)
		min_1=np.min(line_len_pop) # 2nd min
		idx1=list(line_len).index(min_1)

		# get the line pts of 1st & 2nd mins

		line0,line1=line_comb[idx0],line_comb[idx1]
		tri_pt_idx=list(set().union(line0,line1))
		#print tri_pt_idx

		if (len(tri_pt_idx)==3):
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

			#print proj_cord

			middle=np.median(proj_cord)
			md_idx=list(proj_cord).index(middle)

			## assign parity
			mid_f=tri_f[md_idx] # this is the flux of img B (should assign negative parity)

			rfold=(np.sum(dou_f)-2.*mid_f)/np.sum(dou_f)
			rcusp=(np.sum(tri_f)-2.*mid_f)/np.sum(tri_f)

			print rfold,rcusp

			#rcusp_file.write(str(rfold)+'\t'+str(rcusp)+'\n')

			# go over subs


		else:
			print "Merging images are too close!"
			drop=drop+1
			rcusp_file.write('# '+str(i)+' source dropped\n')
			drop_log.write(str(i)+'\n')

	else:
		print "Not a quad!"
		drop=drop+1
		rcusp_file.write('# '+str(i)+' source dropped\n')
		drop_log.write(str(i)+'\n')


	

	#if(len(peaks['x_peak'])!=3 & len(peaks['x_peak'])!=4):
	#	print "Not a quad"


print '# of drop ='+str(drop)



