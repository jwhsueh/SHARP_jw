from astropy import units as u
from astropy.coordinates import SkyCoord

coord = ['07:16:03.576 +47:08:50.154','07:16:03.582 +47:08:50.0','07:16:03.656 +47:08:49.49','07:16:03.692 +47:08:50.615']

img = []

for i in coord:
	img.append(SkyCoord(i, unit= (u.deg, u.deg)))

#print img[0].ra.dms
#print img[1].ra.dms
#print img[2].ra.dms
#print img[3].ra.dms

imgRA = []
imgDec = []

RA0 = img[0].ra.degree
Dec0 = img[0].dec.degree

for i in img:
	imgRA.append((i.ra.degree - RA0)*3600.*10.)
	imgDec.append((i.dec.degree - Dec0)*3600.)

print imgRA
print imgDec
