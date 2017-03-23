import numpy as np

def rcusp_tool(img_x,img_y,img_f,len_cen):

	## get rid of the fourth img

	mask = ~np.in1d(img_f,np.min(img_f))
	img_x,img_y,img_f = img_x[mask],img_y[mask],img_f[mask]

	## find image B (to assign parity)
	delta_x=np.max(img_x)-np.min(img_x)
	delta_y=np.max(img_y)-np.min(img_y)

	if (delta_x>delta_y): # projection axis=x
		proj_cord=img_x
	else:
		proj_cord=img_y

	#print proj_cord

	middle=np.median(proj_cord)
	md_idx=list(proj_cord).index(middle)
	md_mask = np.in1d(np.arange(3),md_idx)
	mid_x,mid_y,mid_f= img_x[md_mask],img_y[md_mask],img_f[md_mask]

	## calculate phi0 (delta phi)
	pt_x,pt_y,pt_f = img_x[~md_mask],img_y[~md_mask],img_f[~md_mask]
	#print Table([pt_x,pt_y])
	vec0 = np.array([pt_x[0]-len_cen[0],pt_y[0]-len_cen[1]])
	vec1 = np.array([pt_x[1]-len_cen[0],pt_y[1]-len_cen[1]])
	phi0 = np.arccos(np.dot(vec0,vec1)/np.sqrt(np.dot(vec0,vec0)*np.dot(vec1,vec1))) # rad
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
		

	else:
		dou_x = np.array([pt_x[1],mid_x])
		dou_y = np.array([pt_y[1],mid_y])
		dou_f = np.array([pt_f[1],mid_f])
		

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

	return rfold[0],rcusp[0], phi0, phi1