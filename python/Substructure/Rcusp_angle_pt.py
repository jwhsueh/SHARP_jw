from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from photutils import find_peaks
from astropy.stats import sigma_clipped_stats
from itertools import combinations
from astropy import convolution
from astropy.table import Table
from astropy.wcs import WCS

### change here ###
subID=208021
proj=3
NN='64'
#NN='128'
## sub flag =1
subflag=0

#cflag = 0

# or here 
#real_size =5.21E-06# rad
#real_size = np.degrees(real_size)
real_size = 0.00042214843# in deg
# degree
real_size = real_size*3600*1000 #mas

print real_size


n_src=100 # number of sources
img_size=256
len_cen = np.array([img_size/2-1,img_size/2-1])

mer_beam = 10 # mas
beam = mer_beam/(real_size/img_size)
print beam



if (subflag==0):
	output_name = '/'+str(subID)+'_p'+str(proj)+'_'+NN+'_Rcusp_s.txt' # no sub
	cmask = 50 # for source mask
else:
	output_name = '/'+str(subID)+'_p'+str(proj)+'sub_'+NN+'_Rcusp.txt' # w/ sub
	cmask = 50


#if (cflag==1):
#	output_name = '/'+str(subID)+'_p'+str(proj)+'_'+NN+'_Rcusp_c.txt' # no sub
#	cmask = 50 # for source mask

imagepath='/Volumes/sting_1/snap99_'+str(subID)
filepath='/Volumes/sting_1/snap99_'+str(subID)
outpath='/Volumes/sting_1/data/Rcusp_c'
magpath='/Volumes/sting_1/data/invmag'

mag_filename = '/particles'+str(subID)+'_p'+str(proj)+'_'+NN+'.invmag.fits'
mag_fits=fits.open(filepath+mag_filename)
mag = 1/mag_fits[0].data
mag_wcs = WCS(mag_fits[0].header)


rcusp_file=open(outpath+output_name,'w')
#rcusp_file.write('# n_src = '+str(n_src)+'\n')
#rcusp_file.write('# Rfold\tRcusp\tphi0\tphi1\n')


drop=0
drop_file='/'+str(subID)+'_p'+str(proj)+'_'+NN+'_drop.txt'
drop_log=open(filepath+drop_file,'w')
#drop_log=open(filepath+drop_file,'a+')

for i in range(n_src):
	print "# "+str(i)+" source:"

	
	if (subflag==0):
		filename='/image_'+str(subID)+'_p'+str(proj)+'_'+NN+'src_'+str(i)+'.fits' # no sub
		img_outname= '/image_'+str(subID)+'_p'+str(proj)+'_'+NN+'src_'+str(i)+'_s' # no sub

	else:
		filename='/image_'+str(subID)+'_p'+str(proj)+'sub_'+NN+'src_'+str(i)+'.fits' # w/ sub
		img_outname= '/image_'+str(subID)+'_p'+str(proj)+'sub_'+NN+'src_'+str(i) # w/ sub
	

	#if (cflag==1):
	#	filename='/image_'+str(subID)+'_p'+str(proj)+'_'+NN+'src_'+str(i)+'.fits' # no sub
	#	img_outname= '/image_'+str(subID)+'_p'+str(proj)+'_'+NN+'src_'+str(i)+'_c' # no sub

	image_fits=fits.open(filepath+filename)
	image_wcs=WCS(image_fits[0].header)
	image=image_fits[0].data
	image_or=image
	#plt.imshow(image_or)
	#plt.show()

	## smoothing
	kernel=convolution.Gaussian2DKernel(beam)
	image=convolution.convolve(image,kernel)

	mean,median,std=sigma_clipped_stats(image_or,sigma=3.0)
	threshold=median+(20*std)
	peaks=find_peaks(image,threshold,box_size=5) # w/ convolve 

	if np.max(image) == 0:
		print "No images!"
		drop=drop+1
		#rcusp_file.write('# '+str(i)+' source dropped\n')
		drop_log.write(str(i)+'\n')

	#c_peaks=find_peaks(image,threshold,box_size=5) # 

	#print peaks

	else:
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
			
			## turn on when cluster pts are too close
			if (min_0+min_1 < img_size/10):
				
				# replace the fainter one in line0 as the fourth img
				pt0,pt1 = line0[0],line0[1]
				fourth_idx = quad_idx[~tri_mask][0]
				if (lens_f[pt0] > lens_f[pt1]):
					fake_idx = pt1
				else:
					fake_idx = pt0

				line0[line0==fake_idx] = fourth_idx
				line1[line1==fake_idx] = fourth_idx
			

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
				## re-calculate line length of line0, line1
				pt0,pt1 = line0[0],line0[1]
				len0 = np.sqrt((lens_x[pt0]-lens_x[pt1])**2+(lens_y[pt0]-lens_y[pt1])**2)
				pt0,pt1 = line1[0],line1[1]
				len1 = np.sqrt((lens_x[pt0]-lens_x[pt1])**2+(lens_y[pt0]-lens_y[pt1])**2)

				if (len1 < len0):
					line_long,line_short = line0,line1
					line0,line1 = line_short,line_long
			
				## find magnification
				
				#for i in range(3):
					#frame_y,frame_x = np.ogrid[-tri_y[i]:img_size-tri_y[i],-tri_x[i]:img_size-tri_x[i]]
					#mask = frame_x*frame_x+frame_y*frame_y < (beam)**2
					#tri_f[i] = np.sum(image[mask])
				peak_wcs = image_wcs.wcs_pix2world(tri_x,tri_y,1)
				px,py = mag_wcs.wcs_world2pix(peak_wcs[0],peak_wcs[1],1)
				print px,py
				for i in range(3):
					tri_f[i] = mag[int(py[i]),int(px[i])]

				tri_table = Table([tri_x,tri_y,tri_f])
				print tri_table
				'''
				## save merging double index
				dou_mask=np.in1d(quad_idx,line0)
				dou_x,dou_y,dou_f=lens_x[dou_mask],lens_y[dou_mask],lens_f[dou_mask]
				
				#for i in range(2):
					#frame_y,frame_x = np.ogrid[-dou_y[i]:img_size-dou_y[i],-dou_x[i]:img_size-dou_x[i]]
					#mask = frame_x*frame_x+frame_y*frame_y < (beam)**2
					#dou_f[i] = np.sum(image[mask])
				peak_wcs = image_wcs.wcs_pix2world(dou_x,dou_y,1)
				px,py = mag_wcs.wcs_world2pix(peak_wcs[0],peak_wcs[1],1)
				for i in range(2):
					dou_f[i] = mag[int(py[i]),int(px[i])]

				## --- opening angle
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
				md_mask = np.in1d(np.arange(3),md_idx)

				## calculate phi0 (delta phi)
				pt_x,pt_y = tri_x[~md_mask],tri_y[~md_mask]
				#print Table([pt_x,pt_y])
				vec0 = np.array([pt_x[0]-len_cen[0],pt_y[0]-len_cen[1]])
				vec1 = np.array([pt_x[1]-len_cen[0],pt_y[1]-len_cen[1]])
				phi0 = np.arccos(np.dot(vec0,vec1)/np.sqrt(np.dot(vec0,vec0)*np.dot(vec1,vec1))) # rad
				phi0 = np.degrees(phi0)

				## calculate phi1
				pt_x,pt_y = dou_x,dou_y
				#print Table([pt_x,pt_y])
				vec0 = np.array([pt_x[0]-len_cen[0],pt_y[0]-len_cen[1]])
				vec1 = np.array([pt_x[1]-len_cen[0],pt_y[1]-len_cen[1]])
				phi1 = np.arccos(np.dot(vec0,vec1)/np.sqrt(np.dot(vec0,vec0)*np.dot(vec1,vec1))) # rad
				phi1 = np.degrees(phi1)		

				## assign parity
				mid_f=tri_f[md_idx] # this is the flux of img B (should assign negative parity)

				rfold=(np.sum(dou_f))/np.sum(np.abs(dou_f))
				rcusp=(np.sum(tri_f))/np.sum(np.abs(tri_f))

				print rfold,rcusp,phi0,phi1
				'''
				#rcusp_file.write(str(rfold)+'\t'+str(rcusp)+'\t'+str(phi0)+'\t'+str(phi1)+'\n')

				# closest triplet img
				plt.imshow(image)
				plt.scatter(tri_x,tri_y,marker='o',s=100,edgecolor='r',facecolor='none')
				#plt.show()
				plt.savefig(imagepath+img_outname+'.png')
				plt.clf()

				# closest triplet txt
				tri_table = Table([tri_x,tri_y,tri_f])
				tri_table.write(imagepath+img_outname+".txt",format='ascii')


			else:
				print "Merging images are too close!"
				drop=drop+1
				#rcusp_file.write('# '+str(i)+' source dropped\n')
				drop_log.write(str(i)+'\n')

		else:
			print "Not a quad!"
			drop=drop+1
			#rcusp_file.write('# '+str(i)+' source dropped\n')
			drop_log.write(str(i)+'\n')


	

	#if(len(peaks['x_peak'])!=3 & len(peaks['x_peak'])!=4):
	#	print "Not a quad"


print '# of drop ='+str(drop)



