import numpy as np

path = '/Volumes/narsil_1/jwhsueh/pylens_mc/'
lens = 'MG0414'
fsub = '0133'

rayfile = open(path+lens+'/'+lens+'_'+fsub+'_'+'1out.txt','r')
lines = rayfile.readlines()
#lines = lines[:50]

outfile = open(path+lens+'/'+lens+'_'+fsub+'_'+'out_com.txt','a')

## first line
fline = lines[0].split()
lines = lines[1:]
quad_list = []
quad_list.append(fline[1:])

for Aline in lines:
	#print quad_list
	Aline = Aline.split()
	Aline[-1] = Aline[-1].split(']')[0]
	#print Aline
	if Aline[0]=='[':

		if len(quad_list)==4:
			## orientation!
			quad_list = np.array(quad_list).astype(float)
			outfile.write(str(quad_list[:,0])[1:-1]+'\t'+str(quad_list[:,1])[1:-1]+'\t'+str(quad_list[:,2])[1:-1]+'\t'+str(quad_list[:,3])[1:-1]+'\n')
			#outfile.write(str(quad_list[0,:3])[1:-1]+'\t'+str(quad_list[:3,1])[1:-1]+'\t'+str(quad_list[:3,2])[1:-1]+'\t'+str(quad_list[:3,3])[1:-1]+'\n')
		quad_list = []
		quad_list.append(Aline[1:])

	else:
		quad_list.append(Aline)
