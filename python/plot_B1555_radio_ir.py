"""
A program to generate figure 2 in the B1555 edge-on disk paper.
This figure shows the Keck AO image as greyscale, and then has contours
that show the MERLIN emission.

Typical output file name: 1555_ao_merlin_overlay.eps
"""

import imfuncs as imf
from matplotlib import pyplot as plt
import sys

""" Set the desired output """
if len(sys.argv)>1:
    outname = sys.argv[1]
else:
    outname = None

""" Set the input file names """
aoim      = 'B1555_nirc2_n_Kp_6x6.fits'
merlin_im = '1555_merlin_5ghz.fits'
vlbi_im   = 'B1555_NA_SELFCAL.FITS'

""" Hardwire rms levels if needed """
rms_vlbi = 0.0001

""" Set the image center and size """
racent  = 239.29968
deccent = 37.359921
imsize  = 1.2       # Value in arcsec

""" Set the figure size """
plt.figure(figsize=(5.7,5.7))

""" Make the overlay plot """
imf.overlay_contours(merlin_im,vlbi_im,racent,deccent,imsize,rms2=rms_vlbi,
                     showradec=False,sighigh=6.)

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

            

""" Label the lensed images """
plt.text(0.33,0.30, 'A',fontdict=font)
plt.text(0.17,0.36, 'B',fontdict=font)
plt.text(-0.25,0.25,'C',fontdict=font)
plt.text(0.16,-0.24,'D',fontdict=font)

""" Show/save the figure """
if outname:
    plt.savefig('%s' % outname)
    print ''
    print 'Saved figure to %s' % outname
    print ''
else:
    plt.show()

