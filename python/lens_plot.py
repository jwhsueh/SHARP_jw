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

"""key in system name, # of lens component, output format"""

# system name
lens_name=raw_input('Lens system:')

# number of component
n_com=raw_input('# of lens component:')
n_com=int(n_com)

# Critical curve plot (default) or Image plot (only has obs+model image position)

plot_type=raw_input('0=Image+Critical Curve, 1=Image only:')
plot_type=int(plot_type)






"""Critical curve input file"""
"""Also put that under /data/lens_info dir"""

if plot_type==0:
    file_name=raw_input('file name of _crit.dat file:')
    name='../models/lens_info/'+file_name+'_crit.dat'



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

##----read findimg positions from files [glafic]----##

def read_findimg(infile):
	## read gravlens image positions
	x,y=np.loadtxt(infile,unpack=True,usecols=(0,1))

	return x,y


##----read lens model infro from gravlens output 'best.dat'-----##

def read_lens_gravlens(lens_name):

    optfile = '../models/lens_info/'+lens_name+'_best.dat'

    best=open(optfile,'r')

    lines=best.read().splitlines()

    # component read in (centroid position)
    cx,cy = np.zeros(n_com), np.zeros(n_com)
    for i in range(n_com):
        paras = lines[i+1].split()

        cx[i],cy[i] = float(paras[2]), float(paras[3])

    # source position

    paras = lines[2+n_com].split()
    print paras

    sx,sy = float(paras[2]), float(paras[3])

    return cx,cy,sx,sy

##----call functions in the end----##

##----

## source & component positions
def model_plot():
    plt.plot(sx,sy,'o',ms=10,mec='k',mfc='r',label='source')
    plt.plot(cx,cy,'^',ms=10,label='lenses',mfc='k')

##----

def img_pos(lens_name):

    obsfile = '../models/lens_info/'+lens_name+'_obs.dat'
    modfile = '../models/lens_info/'+lens_name+'_mod.dat' # now is only for gravlens

    x0,y0 = np.loadtxt(obsfile,unpack=True,usecols=(0,1))
    x,y = np.loadtxt(modfile,unpack=True,usecols=(0,1))

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

"""Lens model info (now is only for gravlens)"""
"""put 'best.dat' that under /data/lens_info dir"""

cx,cy,sx,sy = read_lens_gravlens(lens_name)

if plot_type==0:
    pltlm.plot_critcaust(name,'crit')
    pltlm.plot_critcaust(name,'caust',sls=':')

model_plot()
#x0,y0=x,y
#print x,y
img_pos(lens_name)

#plt.xlim(0.4,-0.8)
#plt.ylim(-0.8,0.4)

plt.xlim(2.0,-0.5)
plt.ylim(-1.0,1.5)

plt.legend(loc=1,scatterpoints=1)

if outfile:
	plt.savefig(outfile,bbox_inches='tight')
else:
	plt.show()

