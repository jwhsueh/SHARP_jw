"""
A program to generate figure 2 in the B1555 edge-on disk paper.
This figure shows the Keck AO image as greyscale, and then has contours
that show the MERLIN emission.

Typical output file name: 1555_ao_merlin_overlay.eps
"""

import sys
from matplotlib import pyplot as plt
import numpy as n
import imfuncs as imf

""" Set the desired output """
if len(sys.argv)>1:
    outname = sys.argv[1]
else:
    outname = None

""" Set the figure size """
plt.figure(figsize=(5.7,5.7))

"""
Plots the NIR AO image in greyscale, and then overlays contours
from two radio data sets: MERLIN and VLBA
"""

""" Set the input file names """
aoim      = '../data/B1555_nirc2_n_Kp_6x6.fits'
merlin_im = '../data/1555_merlin_5ghz.fits'
vlbi_im   = '../data/B1555_vlbi_fix_astrom.fits'

""" Hardwire rms levels if needed """
rms_vlbi = 0.0001

""" Set the image center, origin location, and size """
racent  = 239.29968
deccent = 37.359921
zeropos = (0.2276,0.2194)
imsize  = 1.2       # Value in arcsec

""" Make the overlay plot """
imf.overlay_contours(aoim,vlbi_im,racent,deccent,imsize,rms2=rms_vlbi,
                     showradec=False,sighigh=6.,zeropos=zeropos)

""" Set up the font """
if plt.get_backend() == 'MacOSX':
    font = {'style'  : 'normal',
            'color'  : 'black',
            'weight' : 'bold',
            'size'   : 24,
            }
else:
    print 'Using backend %s' % (plt.get_backend())
    font = {'family' : 'serif',
            'style'  : 'normal',
            'color'  : 'black',
            'weight' : 'bold',
            'size'   : 24,
            }

""" 
Label the lensed images, taking into account a possible shift in origin,
which would be set by the zeropos position
"""
labx = n.array([0.33, 0.17, -0.25,  0.16])
laby = n.array([0.30, 0.36,  0.25, -0.24])
labx -= zeropos[0]
laby -= zeropos[1]
labt = ['A', 'B', 'C', 'D']
for i in range(len(labx)):
    plt.text(labx[i],laby[i],labt[i],fontdict=font)

""" Show/save the figure """
if outname:
    plt.savefig('%s' % outname)
    print ''
    print 'Saved figure to %s' % outname
    print ''
else:
    plt.show()

