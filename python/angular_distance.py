### --- change redshifts here ---- ###
z1 = 0.0
z2 = 0.9
### ------------------------------ ###

c = 299792458.
G = 4.3e-6
from math import pi

OMEGA_M,OMEGA_L,h=0.3,0.7,0.7
w_analytic = False
w = -1.
wpars = None

def angular_diameter_distance(z1,z2=0.):
    if z2<z1:
        z1,z2 = z2,z1
    return comoving_transverse_distance(z1,z2)/(1.+z2)

def comoving_transverse_distance(z1,z2=0.):
    dc = 1e5*comoving_distance(z1,z2)/(c/h)
    ok = 1.-OMEGA_M-OMEGA_L
    if ok>0:
        from math import sinh,sqrt
        dtc = sinh(sqrt(ok)*dc)/sqrt(ok)
    elif ok<0:
        from math import sin,sqrt
        ok *= -1.
        dtc = sin(sqrt(ok)*dc)/sqrt(ok)
    else:
        dtc = dc
    return (c/h)*dtc/1e5

def comoving_distance(z1,z2=0.):
    from scipy import integrate
    if z2<z1:
        z1,z2 = z2,z1

    f = lambda z,m,l,k : (m*(1.+z)**3+k*(1.+z)**2+l)**-0.5
    om = OMEGA_M
    ol = OMEGA_L
    ok = 1.-om-ol
    return (c/h)*integrate.quad(f,z1,z2,(om,ol,ok))[0]/1e5


print angular_diameter_distance(z1,z2)