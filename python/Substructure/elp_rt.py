import numpy as np
import matplotlib.pyplot as plt
import commands
import RcuspTool as rt

model_path = '/Users/jwhsueh/Documents/SHARP_jw/data/sub_gravlens/snap99_elp/'

sublist='285514_p1'
n_src=100

## read-in best model
bestfile = open(sublist+'_best.dat','r')
best_line = []
for line in bestfile.readlines():
	best_line.append(line)

lens_model = best_line[1]
model_para = lens_model.split()
lens_cen = np.array([float(model_para[2]),float(model_para[3])])
#print lens_cen

## output crit curve
crit_file = open('crit.input','w')
crit_file.write('gridmode 1\nset ngrid1=25\nset ngrid2=25\nset maxlev=3\n')
crit_file.write('startup 1 1\n')
crit_file.write(lens_model+'\n')
crit_file.write('0 0 0 0 0 0 0 0 0 0\n')
crit_file.write('plotcrit '+sublist+'_crit.dat')
crit_file.close()

opt_out = commands.getstatusoutput('./lensmodel crit.input')
print opt_out[1]


## read in critical curve

tab = np.loadtxt(sublist+'_crit.dat')
ux,uy = np.append(tab[:,2],tab[:,6]),np.append(tab[:,3],tab[:,7])

# dist to lens center
u_dist = np.sqrt((ux-lens_cen[0])**2+(uy-lens_cen[1])**2)
#u_mid = np.median(u_dist)
u_mid = np.average(u_dist)

# 1st & 3rd quad

u_sm = u_dist[u_dist<=u_mid]
u_la = u_dist[u_dist>u_mid]

#u1,u3 = np.median(u_sm),np.median(u_la)
u1,u3 = np.average(u_sm),np.average(u_la)

# grouping -> tangential caustic
mask = np.abs(u_dist-u1) <= np.abs(u_dist-u3)
ux,uy,ur = ux[mask],uy[mask],u_dist[mask]

#plt.scatter(ux,uy,marker='.')
#plt.show()

x0,x1,y0,y1 = np.min(ux),np.max(ux),np.min(uy),np.max(uy)
r_min = np.min(u_dist)

## select src pos & ray-tracing

rcusp_out = open(model_path+sublist+'_rcusp.txt','w')
src_out = open(model_path+sublist+'_srclist.txt','w')

sx0,sy0 = 0.0,0.0
i = 0
while i<n_src:
	## randomly sample
	xi = np.random.random_sample()*(x1-x0)+x0
	yi = np.random.random_sample()*(y1-y0)+y0

	# distance to caustic
	ri = np.sqrt((xi-ux)**2+(yi-uy)**2)
	rci = np.sqrt((xi-lens_cen[0])**2+(yi-lens_cen[1])**2)
	r_dist = np.min(ri)
	idx = np.argmin(ri)
	rc = ur[idx]

	if np.logical_and(r_dist<0.25*r_min,r_dist>0.01*r_min):
		if (rci<rc):
			if np.logical_and(xi!=sx0,yi!=sy0):
				## do findimg

				findimg_file = open('findimg.input','w')
				findimg_file.write('gridmode 1\nset ngrid1=25\nset ngrid2=25\nset maxlev=3\n')
				findimg_file.write('startup 1 1\n')
				findimg_file.write(lens_model+'\n')
				findimg_file.write('0 0 0 0 0 0 0 0 0 0\n')
				findimg_file.write('findimg '+str(xi)+' '+str(yi)+'\n')
				findimg_file.close()

				findimg_out = commands.getstatusoutput('./lensmodel findimg.input')
				print findimg_out[1]

				lines = findimg_out[1].split('\n')
				last_line = lines[-2]
				quad_imghead = lines[-6]
				
				if(quad_imghead[0] == '#'):
					print quad_imghead
					## ** start to change here **
					raytrace_out = open(model_path+'findimg_out.txt','w')
					raytrace_out.write(str(xi)+'\t'+str(yi)+'\n')
					raytrace_out.write(quad_imghead+'\n')
					raytrace_out.write(lines[-5]+'\n'+lines[-4]+'\n'+lines[-3]+'\n'+lines[-2])
					raytrace_out.close()

					## read-in outputs and generate magnification on grid
					tab = np.loadtxt(model_path+'findimg_out.txt',skiprows=1)
					print tab
					img_x,img_y,img_f = tab[:,0],tab[:,1],tab[:,2]

					rfold,rcusp,phi0,phi1 = rt.rcusp_tool(img_x,img_y,img_f,lens_cen)
					print rfold,rcusp,phi0,phi1 
					#rc[i,j],rf[i,j],p0[i,j],p1[i,j] = rcusp,rfold,phi0,phi1

					rcusp_out.write(str(rfold)+'\t'+str(rcusp)+'\t'+str(phi0)+'\t'+str(phi1)+'\n')
					src_out.write(str(xi)+'\t'+str(yi)+'\n')

					sx0,sy0=xi,yi

					i = i+1

				else:
					#rcusp_out.write('# '+str(i)+' drop \n')
					print 'Not a quad'


rcusp_out.close()
src_out.close()
