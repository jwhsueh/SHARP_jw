import pyfits
import numpy as np
import matplotlib.pyplot as plt
from photutils import find_peaks
from astropy.stats import sigma_clipped_stats
from itertools import combinations
from astropy import convolution
from astropy.table import Table

### change here ###
subID=179899
proj=1
NN='64'

rmax =3.69361e-06
real_size = 1.1*2*rmax # rad
# or here
#real_size =1.04336679208202E-05
real_size = np.degrees(real_size)*3600*1000 # mas
print real_size
mer_beam = 10 # mas

n_src=50 # number of sources
img_size=256

beam = mer_beam/(real_size/img_size)
print beam

## sub flag =1
subflag=1

if (subflag==0):
	output_name = '/'+str(subID)+'_p'+str(proj)+'_'+NN+'_Rcusp.txt' # no sub
	cmask = 30 # for source mask
else:
	output_name = '/'+str(subID)+'_p'+str(proj)+'sub_'+NN+'_Rcusp.txt' # w/ sub
	cmask = 50

imagepath='/Volumes/sting_1/snap99_'+str(subID)
filepath='/Volumes/sting_1/snap99_'+str(subID)
outpath='/Volumes/sting_1/data/Rcusp'


rcusp_file=open(outpath+output_name,'w')
rcusp_file.write('# n_src = '+str(n_src)+'\n')
rcusp_file.write('# Rfold\tRcusp\n')


drop=0
drop_file='/'+str(subID)+'_p'+str(proj)+'_'+NN+'_drop.txt'
#drop_log=open(filepath+drop_file,'w')
drop_log=open(filepath+drop_file,'a+')

for i in range(n_src):
	print "# "+str(i)+" source:"

	if (subflag==0):
		filename='/image_'+str(subID)+'_p'+str(proj)+'_'+NN+'src_'+str(i)+'.fits' # no sub
		img_outname= '/image_'+str(subID)+'_p'+str(proj)+'_'+NN+'src_'+str(i)+'.png' # no sub
	else:
		filename='/image_'+str(subID)+'_p'+str(proj)+'sub_'+NN+'src_'+str(i)+'.fits' # w/ sub
		img_outname= '/image_'+str(subID)+'_p'+str(proj)+'sub_'+NN+'src_'+str(i)+'.png' # w/ sub

	image_fits=pyfits.open(filepath+filename)
	image=image_fits[0].data
	image_or=image

	## smoothing
	kernel=convolution.Gaussian2DKernel(beam)
	image=convolution.convolve(image,kernel)

	mean,median,std=sigma_clipped_stats(image,sigma=3.0)
	threshold=median+(20*std)
	peaks=find_peaks(image,threshold,box_size=5)

	## ---- In case that source sometimes show up on the grid

	dist2cen = np.sqrt((peaks['x_peak']-img_size/2)**2+(peaks['y_peak']-img_size/2)**2)
	mask = list(dist2cen).index(min(dist2cen))
	#print min(dist2cen)
	if (dist2cen[mask]<cmask):
		print "mask out the source (central part)"
		peaks.remove_row(mask)

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

		peaks.add_row([corner_x[idx],corner_y[idx],-1])
		lens_x,lens_y,lens_f=peaks['x_peak'],peaks['y_peak'],peaks['peak_value']

	'''
	if (len(peaks['x_peak']) == 5):
		print "mask out the source"
		mask = list(peaks['peak_value']).index(min(peaks['peak_value']))
		peaks.remove_row(mask)
		lens_x,lens_y,lens_f=peaks['x_peak'],peaks['y_peak'],peaks['peak_value']
	'''

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
		#	Will reject the combination with too small separation
		## ----------------------

		img_idx=np.arange(len(peaks['x_peak']))
		comb_list=list(combinations(img_idx,4))
		#print comb_list
		comb_dist_sum=np.zeros(len(comb_list))	# np.inf: this combination is abandoned 
		

		for k in range(len(comb_list)):
			active_idx=comb_list[k]
			act_mask=np.in1d(img_idx,active_idx)
			rest_idx=img_idx[~act_mask]

			#print active_idx,rest_idx
			act_x,act_y,act_f=peaks['x_peak'][act_mask],peaks['y_peak'][act_mask],peaks['peak_value'][act_mask]
			rest_x,rest_y,rest_f=peaks['x_peak'][~act_mask],peaks['y_peak'][~act_mask],peaks['peak_value'][~act_mask]

			##  start to reject the ones with too small separation in the combination
			act_pair = list(combinations(active_idx,2))
			#print act_pair
			act_dist = np.zeros(len(act_pair))

			for j in range(len(act_pair)):
				idx0,idx1 = act_pair[j][0], act_pair[j][1]
				act_dist[j] = np.sqrt((peaks['x_peak'][idx0]-peaks['x_peak'][idx1])**2+(peaks['y_peak'][idx0]-peaks['y_peak'][idx1])**2)
			#print np.min(act_dist)

			if (np.min(act_dist)<7):
				comb_dist_sum[k]=np.inf # this set of active pts is abandoned

			## start to find brighter points than each rest points
			# compare maximum peak in active and rest pts pool			
			if(np.max(rest_f)>np.max(act_f)): # rest pts actually contain one of lens image [This is not enough!]
				comb_dist_sum[k]=np.inf
			# this set of active pts is abandoned

			if (comb_dist_sum[k]<np.inf):
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

		#line0,line1=list(line_comb[idx0]),list(line_comb[idx1])
		line0,line1=np.array(line_comb[idx0]),np.array(line_comb[idx1])
		tri_pt_idx=list(set().union(line0,line1)) # here we have the triplet combination
		#print tri_pt_idx
		if (len(tri_pt_idx)>3):
			pt0,pt1 = line0[0],line0[1]
			if (lens_f[pt0] > lens_f[pt1]):
				fake_idx = pt1
			else:
				fake_idx = pt0

			tri_pt_idx.remove(fake_idx)


		tri_mask=np.in1d(quad_idx,tri_pt_idx)
		#print quad_idx, ~tri_mask

		if (min_0+min_1 < img_size/4):
			
			# replace the fainter one in line0 as the fourth img
			pt0,pt1 = line0[0],line0[1]
			fourth_idx = quad_idx[~tri_mask][0]
			if (lens_f[pt0] > lens_f[pt1]):
				fake_idx = pt1
			else:
				fake_idx = pt0

			#print line0,line1
			line0[line0==fake_idx] = fourth_idx
			line1[line1==fake_idx] = fourth_idx
			#print line0,line1

		tri_pt_idx=list(set().union(line0,line1)) # here we have the triplet combination
		#print tri_pt_idx			
		tri_mask=np.in1d(quad_idx,tri_pt_idx)
		tri_x,tri_y,tri_f=lens_x[tri_mask],lens_y[tri_mask],lens_f[tri_mask]

		if (np.in1d(0,tri_x) | np.in1d(img_size-1,tri_x)):
			print "Merging images are too close!"
			drop=drop+1
			rcusp_file.write('# '+str(i)+' source dropped\n')
			drop_log.write(str(i)+'\n')

		elif (len(tri_pt_idx)==3):
			#tri_mask=np.in1d(quad_idx,tri_pt_idx)
			#tri_x,tri_y,tri_f=lens_x[tri_mask],lens_y[tri_mask],lens_f[tri_mask]
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

			rcusp_file.write(str(rfold)+'\t'+str(rcusp)+'\n')

			plt.imshow(image)
			plt.scatter(tri_x,tri_y,marker='*',s=100,color='k')
			plt.savefig(imagepath+img_outname)
			plt.clf()


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



