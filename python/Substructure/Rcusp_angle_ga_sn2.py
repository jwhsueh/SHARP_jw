from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from astropy import convolution
from astropy.modeling import models, fitting 
from astropy.table import Table

listpath='/Volumes/sting_1/data/'
outpath='/Volumes/sting_1/data/shoot_noise'

list_file = listpath+'/snap99_fa_Rcusp.txt'
size_file = listpath+'/snap99_fa_size.txt'

list_tab = np.genfromtxt(list_file,dtype='str')
size_tab = np.loadtxt(size_file)
size,n_src = size_tab[:,0],size_tab[:,1].astype(int)

NN='128'
#NN='64'
mer_beam = 50 # mas
boxsize = 24
img_size=256
len_cen = np.array([img_size/2-1,img_size/2-1])

list_tab = ['179899sie_p1']
size = np.array([0.000150015*2])
#size = np.array([5.49E-04,4.71E-04,4.97E-04,4.72E-04]) # in deg
n_src = np.empty(len(list_tab))
n_src.fill(200)
n_src = n_src.astype(int)

for i in range(len(list_tab)):
	obj_name = list_tab[i]
	list_str = obj_name.split('_p')
	subID, proj = list_str[0],list_str[1]

	filepath='/Volumes/sting_1/data/shoot_noise/'+str(subID)
	#filepath='/Volumes/sting_1/snap99_179899'

	real_size = size[i] # degree
	real_size = real_size*3600*1000 #mas
	beam = mer_beam/(real_size/img_size)
	print 'beam = '+str(beam)

	# read-in drop file
	drop_file = filepath+'/'+obj_name+'_'+NN+'_drop.txt'
	drop_list = np.loadtxt(drop_file)

	output_name = '/'+obj_name+'_'+NN+'_Rcusp_ga.txt' # no sub
	rcusp_file=open(outpath+output_name,'w')
	rcusp_file.write('# n_src = '+str(n_src[i])+'\n')
	rcusp_file.write('# Rfold\tRcusp\tphi0\tphi1\n')

	#n_src[i] = 2

	for j in range(n_src[i]):
		if (~np.in1d(j,drop_list)): # if the src is not dropped

			print '##### \n'+ obj_name+' src '+str(j)+' '+NN+'\n#####'

			# read-in img center
			img_file = filepath+'/image_'+obj_name+'_'+NN+'src_'+str(j)+'_s.txt'
			fits_file='/image_'+obj_name+'_'+NN+'src_'+str(j)+'.fits'
			img_info = np.loadtxt(img_file,skiprows=1)
			img_x,img_y = img_info[:,0],img_info[:,1]
			img_fk,img_f = np.empty(img_x.size),np.empty(img_x.size)


			image_fits=fits.open(filepath+fits_file)
			image=image_fits[0].data
			image_or=image

			## smoothing
			kernel=convolution.Gaussian2DKernel(beam)
			image=convolution.convolve(image,kernel)
			flux_mid = np.median(image)

			## ------ Gaussian fit	
			## img cut
			gaus2d = models.Gaussian2D(amplitude=1., x_mean=0, y_mean=0., x_stddev=1.,y_stddev=1.)
			fit_g = fitting.LevMarLSQFitter()
			
			grid_y,grid_x = np.mgrid[:boxsize,:boxsize]
			grid_yi,grid_xi = np.mgrid[:img_size,:img_size]
			grid_yy,grid_xx = np.mgrid[-1*boxsize:1*boxsize,-1*boxsize:1*boxsize]

			for k in range(img_x.size):
				dis = (grid_xi-img_x[k])**2+(grid_yi-img_y[k])**2
				mask = dis>boxsize**2

				img_cut = np.copy(image)
				img_cut[mask] = 0

				#plt.imshow(img_cut)
				#plt.show()
				#plt.imshow(image)

				
				#img_f[k] = np.sum(img_cut)
				
				img_cut2 = img_cut[img_y[k]-boxsize/2:img_y[k]+boxsize/2,img_x[k]-boxsize/2:img_x[k]+boxsize/2]
				fit_out = fit_g(gaus2d,grid_x,grid_y,img_cut2)
				#max_v = np.max(fit_out(grid_x,grid_y))
				fit_mod = fit_out(grid_x,grid_y)
				#print fit_out.x_mean, fit_out.y_mean
				#print fit_out.x_mean<1
				boxsize_ad = boxsize
				while (fit_out.x_mean < 1):
					boxsize_ad = boxsize_ad-2
					grid_yad,grid_xad = np.mgrid[:boxsize_ad,:boxsize_ad]
					img_cut2 = img_cut[img_y[k]-boxsize_ad/2:img_y[k]+boxsize_ad/2,img_x[k]-boxsize_ad/2:img_x[k]+boxsize_ad/2]
					fit_out = fit_g(gaus2d,grid_xad,grid_yad,img_cut2)
					print boxsize_ad

					if (boxsize_ad<=4):
						break
				
				#plt.imshow(fit_out(grid_xx,grid_yy))
				#plt.title('src '+str(j))
				#plt.colorbar()
				#plt.show()
				#plt.clf()
				
				ga_array = fit_out(grid_xx,grid_yy)
				#image[img_y[k]-boxsize:img_y[k]+boxsize,img_x[k]-boxsize:img_x[k]+boxsize] = image[img_y[k]-boxsize:img_y[k]+boxsize,img_x[k]-boxsize:img_x[k]+boxsize] - ga_array
				#plt.imshow(image)
				#plt.show()
				#plt.imshow(ga_array)
				#plt.show()
				img_fk[k] = fit_out.amplitude.value
				hm = ga_array[ga_array>np.max(ga_array)/2]
				img_f[k] = np.sum(hm)
				

			## find image B (to assign parity)
			delta_x=np.max(img_x)-np.min(img_x)
			delta_y=np.max(img_y)-np.min(img_y)

			if (delta_x>delta_y): # projection axis=x
				proj_cord=img_x
			else:
				proj_cord=img_y

			#print proj_cord

			middle=np.median(proj_cord)
			#print middle
			md_idx=list(proj_cord).index(middle)
			md_mask = np.in1d(np.arange(3),md_idx)
			mid_x,mid_y,mid_f,mid_fk = img_x[md_mask],img_y[md_mask],img_f[md_mask],img_fk[md_mask]
			#print mid_x,mid_y

			## calculate phi0 (delta phi)
			pt_x,pt_y,pt_f,pt_fk = img_x[~md_mask],img_y[~md_mask],img_f[~md_mask],img_fk[~md_mask]
			#print Table([pt_x,pt_y])
			vec0 = np.array([pt_x[0]-len_cen[0],pt_y[0]-len_cen[1]])
			vec1 = np.array([pt_x[1]-len_cen[0],pt_y[1]-len_cen[1]])
			phi0 = np.arccos(np.dot(vec0,vec1)/np.sqrt(np.dot(vec0,vec0)*np.dot(vec1,vec1))) # rad
			#print phi0
			phi0 = np.degrees(phi0)

			## save merging double index
			## -- here
			dou_x,dou_y = np.empty(2),np.empty(2)

			line0 = (pt_x[0]-mid_x)**2+(pt_y[0]-mid_y)**2
			line1 = (pt_x[1]-mid_x)**2+(pt_y[1]-mid_y)**2

			if (line0<line1):
				dou_x = np.array([pt_x[0],mid_x])
				dou_y = np.array([pt_y[0],mid_y])
				dou_f = np.array([pt_f[0],mid_f])
				dou_fk = np.array([pt_fk[0],mid_fk])

			else:
				dou_x = np.array([pt_x[1],mid_x])
				dou_y = np.array([pt_y[1],mid_y])
				dou_f = np.array([pt_f[1],mid_f])
				dou_fk = np.array([pt_fk[1],mid_fk])


			## calculate phi1
			pt_x,pt_y = dou_x,dou_y
			#print Table([pt_x,pt_y])
			vec0 = np.array([pt_x[0]-len_cen[0],pt_y[0]-len_cen[1]])
			vec1 = np.array([pt_x[1]-len_cen[0],pt_y[1]-len_cen[1]])
			phi1 = np.arccos(np.dot(vec0,vec1)/np.sqrt(np.dot(vec0,vec0)*np.dot(vec1,vec1))) # rad
			phi1 = np.degrees(phi1)

			## assign parity
			
			rfold=(np.sum(dou_f)-2.*mid_f)/np.sum(dou_f)
			rcusp=(np.sum(img_f)-2.*mid_f)/np.sum(img_f)

			## cusp correction
			if (phi0<70):
				rcusp = rcusp/2


			print rfold[0],rcusp[0],phi0,phi1

			rcusp_file.write(str(rfold[0])+'\t'+str(rcusp[0])+'\t'+str(phi0)+'\t'+str(phi1)+'\n')

			# closest triplet img
			plt.imshow(image)
			plt.scatter(img_x,img_y,marker='o',s=100,edgecolor='k',facecolor='none')
			#plt.show()
			plt.savefig(filepath+'/image_'+obj_name+'_'+NN+'src_'+str(j)+'_ga.png')
			plt.clf()

			# closest triplet txt
			tri_table = Table([img_x.astype(int),img_y.astype(int),img_f,img_fk],names=['x','y','tot_f','ga_peak'])
			print tri_table
			tri_table.write(filepath+'/image_'+obj_name+'_'+NN+'src_'+str(j)+'_ga.txt',format='ascii')

	rcusp_file.close()
	print '#####'
	print obj_name
	print '#####'




