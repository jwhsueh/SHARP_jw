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
from matplotlib import patches
import plot_lensmod as pltlm

#---------------------------------------------------------------------------

def model_plot(sx, sy, cx, cy):
    plt.scatter(sx,sy,marker='o',s= 50,edgecolor='k',facecolors='r',label='Source')
    #plt.scatter(cx,cy,'^',ms=10,label='lenses',mfc='k')
    plt.scatter(cx[0],cy[0],marker='^',edgecolor='k',s= 100,facecolors='k',label = 'SIE centre')
    plt.scatter(cx[1],cy[1],marker='v',edgecolor='k',s= 100,facecolors='r',label = 'Disc centre')

#---------------------------------------------------------------------------

def img_pos(obsfile, xmod, ymod):

    x0,y0 = n.loadtxt(obsfile,unpack=True,usecols=(0,1))

    plt.scatter(x0,y0,marker ='+',color='b',s= 100,label='Observed')
    plt.scatter(xmod,ymod,marker='o',edgecolor='r',s= 100,facecolors='none',label='Predicted')
    #plt.scatter(sig_dan_fa,DMfrac_dan_fa,edgecolor='b',facecolors = 'none',marker='o',label='morphology pick: Face-on')

#---------------------------------------------------------------------------

def disc_plane(x_sec,y_sec):

    plt.plot(x_sec,y_sec,'k--')
    
#---------------------------------------------------------------------------

def radio_overlay_b0712():
    """
    Plots the NIR AO image in greyscale, and then overlays contours
    from two radio data sets: MERLIN and VLBA
    """

    """ Set the input file names """
    aoim      = '../data/B0712_nirc2_n_Kp_6x6.fits'
    #merlin_im = '../data/B0712_BS251A1.fits'
    vlbi_im   = '../data/B0712_BS251A1.fits'

    """ Hardwire rms levels if needed """

    rms_vlbi = 0.00005

    """ Set the image center, origin location, and size """
    racent  = 109.01522
    deccent = 47.147305
    zeropos = (-0.79375,-0.16978)
    imsize  = 3       # Value in arcsec

    """ Make the overlay plot """
    imf.overlay_contours(aoim,vlbi_im,racent,deccent,imsize,
                         showradec=False,fmax=6.,zeropos=zeropos,rms2=rms_vlbi)
                         #infile3=vlbi_im,rms3=rms_vlbi,ccolor3='b')

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
    labx = n.array([-0.9, -0.8, 0.15,  0.6])
    laby = n.array([-0.20,-0.50,  -1.05, 0.4])
    labx -= zeropos[0]
    laby -= zeropos[1]
    labt = ['A', 'B', 'C', 'D']
    for i in range(len(labx)):
        plt.text(labx[i],laby[i],labt[i],fontdict=font)

#---------------------------------------------------------------------------

def mark_radio_b0712(color='r', radius=0.05, lw=2):
    """
    Plots the NIR AO image in greyscale, and then marks with circles
    the location of the lensed AGN images
    """

    """ Set the input file names """
    aofile = '../data/B0712_nirc2_n_Kp_6x6.fits'

    """ Set the image center, origin location, and size """
    racent  = 109.01522
    deccent = 47.147305
    zeropos = (-0.79375,-0.16978)
    imsize  = 3       # Value in arcsec

    """ Plot the AO image """
    aoim = imf.Image(aofile)
    aoim.display(cmap='gray_inv',subimdef='radec',subimcent=(racent,deccent),
                 subimsize=(imsize,imsize),dispunits='radec',fmax=6.,
                 zeropos=zeropos)

    """ Add circles for radio positions """
    obsfile  = '../models/lens_info/B0712_obs.dat'
    x0,y0 = n.loadtxt(obsfile,unpack=True,usecols=(0,1))
    ax = plt.gca()
    for i in range(x0.size):
        circ = patches.Circle((x0[i],y0[i]),radius,edgecolor=color,
                              facecolor='none',lw=lw)
        ax.add_patch(circ)
    plt.ylabel('')

#---------------------------------------------------------------------------

def gravlens_b0712(ax=None, showylab=True):
    """
    Plots the gravlens model.  This includes the following information:
       observed image positions
       model-predicted image positions
       lens centroids
       critical and caustic curves

    Inputs:
     ax       - the matplotlib axis to which to add the plotted curve.
                The default (None) should be fine for most applications.
    """

    """ Set up input file information """
    #critfile = '../models/B1555/gravlens_code/B1555_expdisk_try5_1.crit'
    critfile = '../models/lens_info/B0712_LO_crit.dat'
    obsfile  = '../models/lens_info/B0712_obs.dat'

    """ Set model-predicted positions """
    xmod = [8.119778e-01 , 1.236877e-04 ,5.593499e-02,1.173966e+00]
    ymod = [-6.630301e-01 ,7.317180e-05,-1.560910e-01 ,4.590439e-01]

    """ Set lens mass centroids """
    cx = [ 7.852072e-01, 8.956939e-01] # x position
    cy = [ 1.424232e-01, 2.004992e-01] # y position

    """ Set the source positions """
    sx,sy= 6.953240e-01,  1.520697e-02

    """ Set the observed disc mid plane """

    #x_sec = [-0.0608,-0.2724]
    #y_sec = [0.3806,-0.8194]

    """ Do the plotting """
    pltlm.plot_critcaust(critfile,'crit',ax=ax)
    pltlm.plot_critcaust(critfile,'caust',sls=':',ax=ax)
    model_plot(sx,sy,cx,cy)
    img_pos(obsfile,xmod,ymod)
    #disc_plane(x_sec,y_sec)

    plt.xlabel(r'$\Delta  \alpha $ (arcsec)')
    if showylab:
        plt.ylabel('$\Delta \delta$ (arcsec)')

    #plt.xlim(0.3724,-0.8276)
    #plt.ylim(-0.8194,0.3806)
    plt.legend(loc=3,scatterpoints=1)
    #plt.axes().set_aspect('equal')

#---------------------------------------------------------------------------

def plot_2panel_b0712():
    """
    Creates a 2-panel figure, with one panel showing the radio/NIR overlay
    and the other showing the results of the lens modeling.
    """

    """ Set x and y limits for the plots """
    x1 =  1.4
    x2 = -0.5
    y1 = -0.7
    y2 =  0.5

    """ Set up the figure to have the correct dimensions """
    plt.figure(figsize=(11.,5.4))

    """ Get rid of the space between the subplots"""
    plt.subplots_adjust(wspace=0.001)

    """ Make the radio overlay plot """
    ax1 = plt.subplot(121)
    radio_overlay_b0712()

    """ Make the lens modeling plot """
    ax2 = plt.subplot(122, sharey=ax1)
    gravlens_b0712(ax=ax2,showylab=False)

    """ Finalize the plotted ranges """
    ax1.set_xlim(x1,x2)
    ax1.set_ylim(y1,y2)
    ax2.set_xlim(x1,x2)
    ax2.set_ylim(y1,y2)

    """ Handle the axis labels"""
    yticklabels = ax2.get_yticklabels()
    #yticklabels = ax1.get_yticklabels() + ax2.get_yticklabels()
    plt.setp(yticklabels, visible=False)

#---------------------------------------------------------------------------

def plot_3panel_b0712():
    """
    Creates a 3-panel figure, with one panel showing the radio/NIR overlay,
    one showing the NIR image with just circles indicating the lensed image
    locations, and the third showing the results of the lens modeling.
    """

    """ Set x and y limits for the plots """
    x1 =  2.0
    x2 = -0.5
    y1 = -1.3
    y2 =  1.2

    """ Set up the figure to have the correct dimensions """
    panelsize = 6.
    xsize = 3.*panelsize + 0.4
    plt.figure(figsize=(xsize,panelsize))

    """ Get rid of the space between the subplots"""
    plt.subplots_adjust(wspace=0.001)

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
                'size'   : 28,
                }

    """ Set up for plot labels """
    lx = -0.1
    ly = -1.2

    """ Make the radio overlay plot """
    ax1 = plt.subplot(131)
    radio_overlay_b0712()
    plt.text(lx,ly,'(a)',fontdict=font)

    """ Make the middle plot """
    ax2 = plt.subplot(132, sharey=ax1)
    mark_radio_b0712()
    plt.text(lx,ly,'(b)',fontdict=font)

    """ Make the lens modeling plot """
    ax3 = plt.subplot(133, sharey=ax1)
    gravlens_b0712(ax=ax3,showylab=False)
    plt.text(lx,ly,'(c)',fontdict=font)

    """ Finalize the plotted ranges """
    ax1.set_xlim(x1,x2)
    ax1.set_ylim(y1,y2)
    ax2.set_xlim(x1,x2)
    ax2.set_ylim(y1,y2)
    ax3.set_xlim(x1,x2)
    ax3.set_ylim(y1,y2)

    """ Handle the axis labels"""
    yticklabels = ax2.get_yticklabels() + ax3.get_yticklabels()
    plt.setp(yticklabels, visible=False)