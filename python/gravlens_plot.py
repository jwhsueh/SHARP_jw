import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpat
import sys
import plot_lensmod as pltlm

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


## inputs for plotting

# input files

name='../models/B1555/gravlens_code/B1555_expdisk_try5_1.crit'
obsfile = '../models/B1555/B1555_obsdat.txt'

#name_findimg='/Users/jwhsueh/Documents/SHARP_jw/DK02/gravlens_mock_findimg.dat'

# glafic

#name_opt='exp_optresult.dat'

# number of components
ncom=2

## NOTE!! I've made change only on x axis of B1555

global x,y,x0,y0,sx,sy,cx,cy

# lensed images

## system == B1555
# model
x=[-0.4117 , -0.1629 ,-0.00,-0.0727]
y=[-0.0281,-0.3653,-0.00 ,0.0477]

# model params
# source
sx,sy= -1.952576e-01, -1.493771e-01

# mass profile centroid
cx=[ -1.818271e-01, -1.471605e-01] # x position
cy=[ -1.987580e-01, -2.056755e-01] # y position
'''


## system == B0712
# observed data
x0=[0.0, -0.075,-1.185,-1.71]
y0=[0.0, -0.16, -0.67,0.46]

# model
x=[-9.694099e-01, 4.912447e-02,-3.756152e-02,-1.223079e+00]
y=[-7.651540e-01,  1.666924e-01, -2.644472e-01, 5.127648e-01]

# model params
# source
sx,sy= -8.809485e-01,  2.514535e-01

# mass profile centroid
cx=[ -8.7e-01, -1.08] # x position
cy=[ 0.15, 0.08] # y position
'''

##----read model details from files [glafic]----##

def read_lens(infile):
	## read glafic lens optimal result
	t=open(infile,'r')

        optfile=[]

        for line in t.readlines():
            optfile.append(line)

        src=optfile[-2]
        comp1=optfile[-4]
        comp2=optfile[-3]

        src=src.split()
        comp1=comp1.split()
        comp2=comp2.split()
    
        sx=float(src[2])
        sy=float(src[3])
        cx=[float(comp1[3]),float(comp2[3])]
        cy=[float(comp1[4]),float(comp2[4])]

	return sx,sy,cx,cy

##----read findimg positions from files [gravlens]----##

def read_findimg(infile):
	## read gravlens image positions
	t=np.loadtxt(infile)

    	x=t[:,0]
    	y=t[:,1]

	return x,y


##----call functions in the end----##

##----

## source & component positions
def model_plot():
    plt.plot(sx,sy,'o',ms=10,mec='k',mfc='r',label='source')
    plt.plot(cx,cy,'^',ms=10,label='lenses',mfc='k')

##----

def img_pos(obsfile):

    x0,y0 = np.loadtxt(obsfile,unpack=True,usecols=(0,1))

    plt.plot(x0,y0,'b+',ms=10,label='observed')
    plt.plot(x,y,'o',ms=10,mec='r',mfc='none',label='predicted')
    
    plt.xlabel(r'$\Delta  \alpha $ (arcsec)')
    plt.ylabel('$\Delta \delta$ (arcsec)')
    #plt.xlabel('arcsec')
    #plt.ylabel('arcsec')




##---- call functions from here ----##

#sx,sy,cx,cy=read_lens(name_opt)
#x,y=read_findimg(name_findimg)

plt.figure(figsize=(5.7,5.7))

pltlm.plot_critcaust(name,'crit')
pltlm.plot_critcaust(name,'caust',sls=':')
model_plot()
#x0,y0=x,y
#print x,y
img_pos(obsfile)

plt.xlim(0.4,-0.8)
plt.ylim(-0.8,0.4)
plt.legend(loc=1)

if outfile:
	plt.savefig(outfile,bbox_inches='tight')
else:
	plt.show()

