"""
A bit of code to create a 3-panel plot that combines the old Figure 2
(radio contours overlaid on AO greyscale) and Figure 3 (lensing caustics
and critical curves) from the submitted version of the B1555 paper into
one figure.
"""

import sys
from matplotlib import pyplot as plt
import plot_fns_b1555 as plt1555

""" Check the command line for an optional output file name """
print ''
if len(sys.argv)>1:
	outfile = sys.argv[1]
	print 'Will save output to %s' % outfile
	print ''
else:
	outfile = None
	print 'No output file requested.'
	print ''

""" Create the 3-panel figure """
plt1555.plot_3panel_b1555()

""" Save the output, if requested """
if outfile:
	plt.savefig(outfile,bbox_inches='tight')
else:
	plt.show()
