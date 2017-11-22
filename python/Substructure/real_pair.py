import numpy as np

#pair_list = open('/Volumes/sting_1/B1422_0100_pair.txt','w')

sub_idx = np.arange(110000)
los_idx = np.arange(10000)

n_pair = 110000

sub_pair = np.random.choice(sub_idx,n_pair)
los_pair = np.random.choice(los_idx,n_pair)

np.savetxt('/Volumes/sting_1/subs/B1422_0100_pair.txt', np.c_[sub_pair,los_pair],fmt='%d')