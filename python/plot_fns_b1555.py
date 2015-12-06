"""
A collection of functions that are used to make the plots for B1555.

These are split off from their original python programs to allow for
the possibility of doing a multi-panel plot for the B1555 paper.

Functions
 radio_overlay_b1555  - Overlays contours from MERLIN and VLBA onto AO image
 gravlens_b1555       - Plots gravlens model
"""

import numpy as n
import imfuncs as imf
from matplotlib import pyplot as plt
import plot_lensmod as pltlm

#---------------------------------------------------------------------------
def model_plot(sx, sy, cx, cy):
    plt.plot(sx,sy,'o',ms=10,mec='k',mfc='r',label='source')
    plt.plot(cx,cy,'^',ms=10,label='lenses',mfc='k')

#---------------------------------------------------------------------------

def img_pos(obsfile, xmod, ymod, showylab=True):

    x0,y0 = n.loadtxt(obsfile,unpack=True,usecols=(0,1))

    plt.plot(x0,y0,'b+',ms=10,label='observed')
    plt.plot(xmod,ymod,'o',ms=10,mec='r',mfc='none',label='predicted')
    
    plt.xlabel(r'$\Delta  \alpha $ (arcsec)')
    if showylab:
        plt.ylabel('$\Delta \delta$ (arcsec)')

#---------------------------------------------------------------------------

def radio_overlay_b1555():
    """
    Set up basic parameters of the B1555 system.
    """

    """ Set the input file names """
    aoim      = 'B1555_nirc2_n_Kp_6x6.fits'
    merlin_im = '1555_merlin_5ghz.fits'
    vlbi_im   = 'B1555_vlbi_fix_astrom.fits'

    """ Hardwire rms levels if needed """
    rms_vlbi = 0.0001

    """ Set the image center, origin location, and size """
    racent  = 239.29968
    deccent = 37.359921
    zeropos = (0.2236,0.2174)
    imsize  = 1.2       # Value in arcsec

    """ Make the overlay plot """
    imf.overlay_contours(aoim,merlin_im,racent,deccent,imsize,
                         showradec=False,sighigh=6.,zeropos=zeropos,
                         infile3=vlbi_im,rms3=rms_vlbi,ccolor3='b')

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

#---------------------------------------------------------------------------

def gravlens_b1555():
    """
    Plots the gravlens model.  This includes the following information:
       observed image positions
       model-predicted image positions
       lens centroids
       critical and caustic curves
    """

    """ Set up input file information """
    critfile = '../models/B1555/gravlens_code/B1555_expdisk_try5_1.crit'
    obsfile  = '../models/B1555/B1555_obsdat.txt'

    """ Set model-predicted positions """
    xmod = [-0.4117 , -0.1629 ,-0.00,-0.0727]
    ymod = [-0.0281,-0.3653,-0.00 ,0.0477]

    """ Set lens mass centroids """
    cx = [ -1.818271e-01, -1.471605e-01] # x position
    cy = [ -1.987580e-01, -2.056755e-01] # y position

    """ Set the source positions """
    sx,sy= -1.952576e-01, -1.493771e-01

    """ Do the plotting """
    pltlm.plot_critcaust(critfile,'crit')
    pltlm.plot_critcaust(critfile,'caust',sls=':')
    model_plot(sx,sy,cx,cy)
    img_pos(obsfile,xmod,ymod)

    plt.xlim(0.4,-0.8)
    plt.ylim(-0.8,0.4)
    plt.legend(loc=1)
    plt.axes().set_aspect('equal')

#---------------------------------------------------------------------------

def plot_2panel_b1555():
    """
    Creates a 2-panel figure, with one panel showing the radio/NIR overlay
    and the other showing the results of the lens modeling.
    """
