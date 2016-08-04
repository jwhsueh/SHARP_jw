import os
import numpy as np
import commands
import gravlens_tool as gTool

### actually re-do optimize

path = '../../data/sub_gravlens/'
opt_input = 'rescale_opt.input'

## create opt file

opt_file = open(path+opt_input,'w')
opt_header = open(path+'opt.head','r')

## write head file
head_lines = opt_header.readlines()

for Aline in head_lines:
	opt_file.write(Aline)

## write number of lens and srcs [startup]

# read-in realization files
real_file = 'realization_1.dat'
opt_result = open(path+real_file,'r')

lines = opt_result.readlines()[:-1]
n_lens = len(lines)
findimg_file.write('startup '+str(n_lens)+' 1\n')

## --end of gravlens startup--

## write opt SIE flags

opt_file.write('\n')
opt_file.write('1 1 1 1 1 1 1 0 0 0\n')

## --write a bunch of zeros [pjaffe]--


zeros = '0 0 0 0 0 0 0 0 0 0\n'
for i in range(n_lens-1):
	opt_file.write(zeros)

opt_file.write('\n')

## write opt command

opt_file.write('optimize \n')

opt_file.close()

##------ run opt file

os.system('./lensmodel '+path+opt_input)


##----- create findimg

find_input = 'rescale_find.input'

findimg_file = open(path+output,'w')

gTool.create_bestSub(path)
opt_result = open(path+'bestSub.dat','r')

lines = opt_result.readlines()[:-1]
n_lens = len(lines)
findimg_file.write('startup '+str(n_lens)+' 1\n')

	# write in substructures

for Aline in lines:

	findimg_file.write(Aline)

opt_result.close()

## --end of gravlens startup--

## --write a bunch of zeros--

zeros = '0 0 0 0 0 0 0 0 0 0\n'
for i in range(n_lens):
	findimg_file.write(zeros)

findimg_file.write('\n')
## write findimg command
# get srcs position
opt_result = open('best.dat','r')
raw_line = opt_result.readlines() # get all the rest lines
#print raw_line
raw_line = raw_line[-11].split()
#print raw_line

findimg_file.write('findimg '+raw_line[1]+' '+raw_line[2]+'\n')

opt_result.close()
findimg_file.close()

##------ run findimg

findimg_out = commands.getstatusoutput('./lensmodel '+path+findimg_file)

findimg_out=findimg_out[1].split('\n')

	#print findimg_out

result = open(path+'findimg.out','w')

	## checking img #

check_line = findimg_out[-6].split()
	# if it's four images
if check_line[0] == '#':
		# save x y mag
		#x,y,mag = [],[],[]
	for i in [-2,-3,-4,-5]:
		Aline = findimg_out[i].split()
			#x.append(float(Aline[0]))
			#y.append(float(Aline[1]))
			#mag.append(float(Aline[2]))
		result.write(Aline[0]+' '+Aline[1]+' '+Aline[2]+'\n')

		flag = True

		#result = [x,y,mag] # write into a file
else:
	print '* Not a quad-system'
	flag = False
	result.write('# Not a quad-system')

result.close()

if flag == True:
	