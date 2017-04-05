import numpy as np
import DistanceTool as distance

zs = 1.34
zl = 0.41
rl = 0.64 # arcsec

class cosmopara:
	OM=0.27
	h=0.7

filename = '../models/B0712/B0712_env.txt'

tab = np.loadtxt(filename)
zp = tab[:,1]
log_Ms = tab[:,2]
theta_s = tab[:,0] # angular separation, arcsec

mask = zp>=0
zp,log_Ms,theta_s = zp[mask],log_Ms[mask],theta_s[mask]

# log sigma
log_sig = 0.301*(log_Ms-10.26)+2.29
sig = np.power(10,log_sig) # km/s
print sig

# Einstein radius
dp = np.zeros(len(zp))
for i in range(len(zp)):
	dp = distance.angular_distance(cosmopara,zp[i])

dp = distance.angular_distance(cosmopara,0.2906)
sig = 398
ds = distance.angular_distance(cosmopara,zs)
dl = distance.angular_distance(cosmopara,zl)
dps = ds-dp

rp = 4.0*np.pi*(sig/3e5)**2*dps/ds # radian
rp = np.degrees(rp)*3600 # arcsec

print rp
beta = (dl-dp)*ds/(dl*(ds-dp))
print beta

# beta factor
beta = np.zeros(len(zp))

for i in range(len(zp)):
	if zp[i]>zl:
		dp = distance.angular_distance(cosmopara,zp[i])

		beta[i] = (dp-dl)*ds/(dp*(ds-dl))

# flexion
f = (1.0-beta)**2*(rl*rp)**2/theta_s**3
print f