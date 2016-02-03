import numpy as np
import scipy.integrate
import matplotlib.pyplot as plt

# power law index, Sigma =  Sigma_0 * r**N
N_bulge = -0.87
N_disk = -0.4

# Sigma_0 for bulge & disk component (at r = 1 kpc)
S_bulge = 5.0
S_disk = 2.0

radius = 10.0 # integrate to 10 kpc

'''
# power law checking figure

p =  np.logspace(-2,1,10000)

k_bulge = S_bulge*p**N_bulge
k_disk = S_disk*p**N_disk

plt.loglog(p,k_bulge)
plt.loglog(p,k_disk,'r')
plt.ylim(0.1,1000)
plt.xlim(10**(-1.5),10)

plt.ylabel('$\Sigma$ / $\Sigma_{crit}$')
plt.xlabel('Radius/kpc')

plt.show()


'''
# Abel transform
def Abel_trans(s, N):
	#F=scipy.integrate.quad(lambda r: r**(N+1)/np.sqrt(r**2-s**2), s, 10)
	F = 

	return 2.0*F[0]

def I_circle(N):
	I0 = scipy.integrate.quad(lambda s: Abel_trans(s,N)*s, 0, radius)

	return 2.0*np.pi*I0[0]

K_bulge = I_circle(N_bulge)*S_bulge
K_disk = I_circle(N_disk)*S_disk

f_disk = K_disk/(K_disk+K_bulge)

print f_disk
#'''