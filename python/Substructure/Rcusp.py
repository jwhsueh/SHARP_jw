import pyfits
import numpy as np
import matplotlib.pyplot as plt
from photutils import find_peaks
from astropy.stats import sigma_clipped_stats
from itertools import combinations
from astropy import convolution

filepath='/Users/jwhsueh/Documents/glamer/examples/ParticleExample/build/'
filename='image_281185sub32src_0.fits'

image_fits=pyfits.open(filepath+filename)
image=image_fits[0].data
image_or=image

## smoothing
kernel=convolution.Gaussian2DKernel(5)
image=convolution.convolve(image,kernel)

# get image size
img_size=image.shape[0]
print img_size

mean,median,std=sigma_clipped_stats(image,sigma=3.0)
threshold=median+(20*std)
peaks=find_peaks(image,threshold,box_size=5)

print peaks
'''
## throw out source peak
# by distance to center
center=img_size/2
dist_center=np.sqrt((peaks['x_peak']-center)**2+(peaks['y_peak']-center)**2)
print dist_center
idx=list(dist_center).index(np.min(dist_center))

threshold=np.min(peaks['peak_value'][idx])

peaks=find_peaks(image,threshold,box_size=5)
print peaks
#print len(peaks['x_peak'])
'''
if (len(peaks['x_peak']) <4):
	print "Number of lensed images is less than four!!"

#plt.imshow(image_or)
#plt.scatter(peaks['x_peak'],peaks['y_peak'],marker='*',s=100,color='k')
#plt.show()

'''

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

for i in range(len(comb_list)):
	active_idx=comb_list[i]
	act_mask=np.in1d(img_idx,active_idx)
	rest_idx=img_idx[~act_mask]

	#print active_idx,rest_idx
	act_x,act_y,act_f=peaks['x_peak'][act_mask],peaks['y_peak'][act_mask],peaks['peak_value'][act_mask]
	rest_x,rest_y,rest_f=peaks['x_peak'][~act_mask],peaks['y_peak'][~act_mask],peaks['peak_value'][~act_mask]

	## start to find brighter points than each rest points
	# compare maximum peak in active and rest pts pool
	
	if(np.max(rest_f)>np.max(act_f)): # rest pts actually contain one of lens image [This is not enough!]
		comb_dist_sum[i]=1e10
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

		comb_dist_sum[i]=dist_sum/np.sum(act_f)
	
	#print np.max(rest_f),np.max(act_f)
	#print comb_dist_sum[i]

print comb_dist_sum
print np.min(comb_dist_sum)
sur_comb_idx = list(comb_dist_sum).index(np.min(comb_dist_sum))
sur_mask = np.in1d(img_idx,comb_list[sur_comb_idx])

lens_x,lens_y,lens_f=peaks['x_peak'][sur_mask],peaks['y_peak'][sur_mask],peaks['peak_value'][sur_mask]

print lens_x,lens_y,lens_f
'''

lens_x,lens_y,lens_f=peaks['x_peak'],peaks['y_peak'],peaks['peak_value']

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

plt.imshow(image_or)
plt.scatter(tri_x,tri_y,marker='*',s=100,color='k')
plt.show()

'''
plt.imshow(image)
plt.savefig('../../data/glamer/test0_img.png')
plt.clf()

plt.imshow(image)
plt.scatter(peaks['x_peak'],peaks['y_peak'],marker='*',s=100,color='k')
plt.savefig('../../data/glamer/test0_peaks.png')
plt.clf()

plt.imshow(image)
plt.scatter(lens_x,lens_y,marker='*',s=100,color='k')
plt.savefig('../../data/glamer/test0_quad.png')
plt.clf()

plt.imshow(image)
plt.scatter(tri_x,tri_y,marker='*',s=100,color='k')
plt.savefig('../../data/glamer/test0_cusp.png')
plt.clf()

plt.imshow(image)
plt.scatter(dou_x,dou_y,marker='*',s=100,color='k')
plt.savefig('../../data/glamer/test0_fold.png')
plt.clf()
#plt.show()
'''