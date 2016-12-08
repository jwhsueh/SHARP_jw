"""
A program to generate figure 2 in the B1555 edge-on disk paper.
This figure shows the Keck AO image as greyscale, and then has contours
that show the MERLIN emission.

Typical output file name: 1555_ao_merlin_overlay.eps
"""

import sys
from matplotlib import pyplot as plt
from plot_fns_B1555 import radio_overlay_b1555

""" Set the desired output """
if len(sys.argv)>1:
    outname = sys.argv[1]
else:
    outname = None

""" Set the figure size """
plt.figure(figsize=(5.7,5.7))

""" Make the figure """
radio_overlay_b1555()

""" Show/save the figure """
if outname:
    plt.savefig('%s' % outname)
    print ''
    print 'Saved figure to %s' % outname
    print ''
else:
    plt.show()

