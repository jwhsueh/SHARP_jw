import  threading, datetime, sys, time
import numpy as np 
#import cluster_emcee as mcmc
import realization_rt as rt
from Queue import Queue

lens = 'B1422'
path_sub = './'+lens+'/real_01/'
path_out = './'+lens+'/'

N_sub = 5000 # number of sub realization
zlens = 0.34
zsrc = 3.62
idx_a = 2 # idx of img A

#def thread_job(idx_sub,idx_th,l_array):
#	likelihood_th = rt.gravlens_rt(idx_sub,path_sub,path_out,lens,zlens,zsrc,idx_th)

#	l_array

#	return likelihood_th

# Objects ==============================================================  
class Job:    
     def __init__(self, idx_sub):    
        self.idx = idx_sub    
     def do(self):    
		#mcmc.emcee_sub(self.idx,path_sub,path_out,lens,zlens,zsrc,idx_a)
        rt.gravlens_rt(self.idx,path_sub,path_out,lens,zlens,zsrc,idx_a)
		#np.savetxt(path_sub+'likelihood_real'+str(self.idx),like_real)

def doJob(*args):  
     queue = args[0]  
     print('thread running')
     while queue.qsize() > 0:  
          job = queue.get()  
          job.do()  


# Paramaters ===========================================================
tot_thread=int(sys.argv[1])
tot_job=int(sys.argv[2])
job_start = int(sys.argv[3])

          
# main =================================================================  
# set up queue  
que = Queue()  
for i in range(tot_job):  
        que.put(Job(i+job_start))  
	print("\t[Info] Queue size={0}...".format(que.qsize())) 
  
# Open threads  
thd=[]
for n in range(0,tot_thread):
    thd.append(threading.Thread(target=doJob, name='Thd'+str(n), args=(que,))) 
  
# Start activity to digest queue.  
st = datetime.datetime.now()  
for n in range(0,tot_thread):
    thd[n].start()  
'''
# Wait for all threads to terminate.  
while True:
    alive=sum([thd_n.is_alive() for thd_n in thd])  
    if alive!=0:
        time.sleep(0.1)
    else: 
        break
# print out running time for check multithreading  
td = datetime.datetime.now() - st  
print("\t[Info] Spending time={0}!".format(td))  
'''