import numpy as np
import angular_distance as DA

### change the file path here ###
file_path = '/Volumes/narsil_1/jwhsueh/illustris_1/snapshot89_re.dat' # path of the input table
output_path = '/Volumes/narsil_1/jwhsueh/illustris_1/' # path of the output table
###---------------------------###

table = np.loadtxt(file_path)
subid = table[:,0]
Re = table[:,1:]

z_illus = np.array([0.9,2.0])
z_b2045 = np.array([0.86,2.35])

### calculate conversoin factor
## D = Dls/Dl/Ds
## Re ~ sqrt(D)
###

D_illus = DA.angular_diameter_distance(z_illus[0],z_illus[1])/DA.angular_diameter_distance(z_illus[0])/DA.angular_diameter_distance(z_illus[1])
D_b2045 = DA.angular_diameter_distance(z_b2045[0],z_b2045[1])/DA.angular_diameter_distance(z_b2045[0])/DA.angular_diameter_distance(z_b2045[1])

Re = Re*np.sqrt(D_b2045/D_illus)

np.savetxt(output_path+'snap89_B2045_re.dat',np.c_[subid,Re],fmt='%d %f %f %f')