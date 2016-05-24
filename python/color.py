import pyfits
import numpy as np
import matplotlib.pyplot as plt

dir_path = '../models/B0712'

Kp_band = pyfits.open(dir_path+'/B0712_nirc2_rescale.fits')
F814W = pyfits.open(dir_path+'/B0712_wfpc2_F814W_rescale.fits')
F555W = pyfits.open(dir_path+'/B0712_wfpc2_F555W_rescale.fits')

Kp_img = Kp_band[0].data
F814W_img = F814W[0].data
F555W_img = F555W[0].data

##### color maps

blue = F555W_img/F814W_img
red = F814W_img/Kp_img

##### create mask

lens_pos = np.array([[76,75,60,54],[57,54,44,67]])
r = 3

img_x = np.arange(Kp_img.shape[0])
img_y = np.arange(Kp_img.shape[1])
x_mesh,y_mesh = np.meshgrid(img_x,img_y)

def distance(pos_index):
	dist_array = np.zeros(Kp_img.shape)

	#print lens_pos[0,pos_index],lens_pos[1,pos_index]
	dist_x = x_mesh - lens_pos[0,pos_index]
	dist_y = y_mesh - lens_pos[1,pos_index]

	dist_array = np.sqrt(dist_x**2+dist_y**2)

	return dist_array

class dist:
	A,B,C,D = distance(0),distance(1),distance(2),distance(3)

class mask:
	A,B,C,D = [dist.A<=r],[dist.B<=r],[dist.C<=r],[dist.D<=r]

class color1:
	A,B,C,D = np.zeros(Kp_img.shape)*np.nan,np.zeros(Kp_img.shape)*np.nan,np.zeros(Kp_img.shape)*np.nan,np.zeros(Kp_img.shape)*np.nan

	A[mask.A] = blue[mask.A]
	B[mask.B] = blue[mask.B]
	C[mask.C] = blue[mask.C]
	D[mask.D] = blue[mask.D]

class color2:
	A,B,C,D = np.zeros(Kp_img.shape)*np.nan,np.zeros(Kp_img.shape)*np.nan,np.zeros(Kp_img.shape)*np.nan,np.zeros(Kp_img.shape)*np.nan

	A[mask.A] = red[mask.A]
	B[mask.B] = red[mask.B]
	C[mask.C] = red[mask.C]
	D[mask.D] = red[mask.D]

#print Kp_img.shape
#print dist.B
plt.scatter(color1.A,color2.A,edgecolor='r',facecolors = 'none',marker='^',label='A')
plt.scatter(color1.B,color2.B,marker='o',facecolors = 'none',edgecolors='b',label='B')
plt.scatter(color1.C,color2.C,edgecolor='g',facecolors = 'none',marker='^',label='C')
plt.scatter(color1.D,color2.D,edgecolors='k',marker='o',facecolors = 'none',label='D')
plt.legend(scatterpoints=1)
plt.xlabel('F555W/F814W')
plt.ylabel('F814W/Kp')
plt.title('B0712 color-color diagram (r=0.1")')
plt.show()

