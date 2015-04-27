"""
A program to generate figure 2 in the B1555 edge-on disk paper.
This figure shows the Keck AO image as greyscale, and then has contours
that show the MERLIN emission.
"""

import imfuncs as imf
from matplotlib import pyplot as plt

""" Make the overlay plot """
imf.overlay_contours('B1555_nirc2_n_Kp_6x6.fits','1555_merlin_5ghz.fits',
                     239.29968,37.359921,1.2,showradec=False,sighigh=6.)

""" Label the lensed images """
plt.text(0.33,0.30,'A',fontsize=24,color='r')
plt.text(0.17,0.36,'B',fontsize=24,color='r')
plt.text(-0.25,0.25,'C',fontsize=24,color='r')
plt.text(0.16,-0.22,'D',fontsize=24,color='r')

""" Show the figure """
plt.show()
