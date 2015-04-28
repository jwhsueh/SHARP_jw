"""
A program to generate figure 2 in the B1555 edge-on disk paper.
This figure shows the Keck AO image as greyscale, and then has contours
that show the MERLIN emission.
"""

import imfuncs as imf
from matplotlib import pyplot as plt

""" Set the desired output """
output = 'eps'
outname = '1555_ao_merlin_overlay'

""" Set the figure size """
plt.figure(figsize=(5.7,5.7))

""" Make the overlay plot """
imf.overlay_contours('B1555_nirc2_n_Kp_6x6.fits','1555_merlin_5ghz.fits',
                     239.29968,37.359921,1.2,showradec=False,sighigh=6.)

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
if output == 'eps':
    plt.savefig('%s.eps' % outname)
    print ''
    print 'Saved figure to %s.eps' % outname
    print ''
else:
    plt.show()

