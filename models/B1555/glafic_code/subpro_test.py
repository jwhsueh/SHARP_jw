import os
import subprocess as sub
import numpy as np

def call_findimage(paras):

	## create input file
	head=open('input.head','r')
	tail=open('input.tail','r')
	inp=open('emcee.input','w')

	# read head of input template
	for line in head.readlines():
		#print line
		inp.write(line)

	head.close()
	inp.write('\n')

	# prepare for lens model parameters
	lens_para=['','','']

	for i in range(5):
		lens_para[0]=lens_para[0]+'%f '%paras[i]    #bulge
		lens_para[1]=lens_para[1]+'%f '%paras[5+i]  #disk

		if i<2:
			lens_para[2]=lens_para[2]+'%f '%paras[10+i]  #source

    # write lens models into input file
	inp.write('lens sie '+lens_para[0]+'0.0 0.0 \n')
	inp.write('lens sie '+lens_para[1]+'0.0 0.0 \n')
	inp.write('point 1.5 '+lens_para[2]+'\n')

    # read & write tail of input template
	for line in tail.readlines():
		inp.write(line)

	tail.close()
	inp.close()

    ##

    ## run glafic
	os.system('glafic emcee.input') #this line works for linux machine
    
    #sub.call('glafic emcee.input')

    # read findimg result
	r=np.loadtxt('emcee_point.dat')

	if len(r[:,0])>1:
		r=r[1:,:] #discard the first line
		r[:,2]=np.absolute(r[:,2]) # ignore parity in magnification
		
		for i in range(len(r[:,0])):
			if i ==0:		
				out=np.array(r[i,:3])
				
			else:
				out=np.append(out,r[i,:3])

	print out
	return out   #return findimg result

		
#########	

    
def lnprob(paras):
    
	model=call_findimage(paras)
    
    # B1555 constraints
	#x= np.arrange(10)
	y= np.array([0.0,0.0,17.0,0.0726,0.0480,10.54,0.4117,-0.0280,8.619,0.1619,-0.3680,1.462])
	err= np.array([1.0e-6,1.0e-6,1.70,0.01,0.01,0.97,0.01,0.01,0.83,0.06,0.06,0.13])

	chi2=0
	for i in range(len(y)):
		chi2=chi2+(model[i]-y[i])**2/err[i]**2

	return -0.5*chi2

###########
p=[150.62905,0.17383945,-0.2428983,0.22057935,91.511635,133.6333,0.15735815,-0.2346221,0.8530964,6.221847,0.17383945,-0.2428983]

#call_findimage(p)
print lnprob(p)
