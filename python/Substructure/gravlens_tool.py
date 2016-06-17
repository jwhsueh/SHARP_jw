import os
import numpy as np
import commands

def create_opt(macro_mod,micro_mod,path):
	opt_file = open(path+'valid_check.input','w')

	## write number of lens and srcs [startup]

	n_lens = 1+len(micro_mod.xi)

	opt_file.write('startup '+str(n_lens)+' 1\n')

	## macro model
	macro_line = gravlens_SIE(macro_mod)
	opt_file.write(macro_line)

	## micro model
	for i in range(len(micro_mod.xi)):
		micro_line = gravlens_pjaffe(micro_mod.bi[i],micro_mod.xi[i],micro_mod.yi[i],micro_mod.rt_i[i])
		opt_file.write(micro_line)

	## --end of gravlens startup--

	## --write a bunch of zeros--

	opt_file.write('\n')
	zeros = '0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0\n'
	for i in range(n_lens):
		opt_file.write(zeros)

	opt_file.write('\n')

	## write opt command

	src_pos = str(macro_mod.src_x)+' '+str(macro_mod.src_y)
	opt_file.write('opt '+src_pos)

	opt_file.close()

	return

def gravlens_SIE(model):

	lenspara = [model.b,model.xc,model.yc,model.q,model.PA,model.gamma1,model.gamma2]
	lenspara = str(lenspara)[1:-1].replace(',','')
	Aline = 'alpha '+lenspara+' 0.0 0.0 1.0\n'

	return Aline

def gravlens_pjaffe(b,x,y,rt):

	lenspara = [b,x,y,0.0,0.0,0.0,0.0,0.0,rt,0.0]
	lenspara = str(lenspara)[1:-1].replace(',','')
	Aline = 'pjaffe '+lenspara[1:-1]+'\n'

	return Aline

def opt_result(opt_file):
	opt_output = commands.getstatusoutput('./lensmodel '+opt_file)

