import numpy as np
import matplotlib.pyplot as plt

dirc='../../data/illustris_1/'

## all
table=np.loadtxt(dirc+'AllGalaxy_099_y.dat',skiprows=1)

subfindID=table[:,0]
df=table[:,1]
bf=table[:,2]
Re=table[:,4]
subflag=table[:,5]
stms=table[:,6]
ser=table[:,8]
mor=table[:,10]

## Re flag
mask=Re>0.
subfindID,bf,stms,ser,mor,subflag,df=subfindID[mask],bf[mask],stms[mask],ser[mask],mor[mask],subflag[mask],df[mask]


##--- central flag
mask=subflag==0.
subfindID,bf,stms,ser,mor,df=subfindID[mask],bf[mask],stms[mask],ser[mask],mor[mask],df[mask]


##

bf_mask=bf<0.6
df_mask=df>0.4
mor_mask=(ser<2.0)&(mor==0.)
tri_mask=(bf<0.6)&(df>0.4)&(ser<2.0)&(mor==0.)

k_id=subfindID[bf_mask]
m_id=subfindID[mor_mask]
tri_id=subfindID[tri_mask]

union=np.union1d(k_id,m_id)
print union[:10]
inter=np.intersect1d(k_id,m_id)

only_k=k_id[~np.in1d(k_id,inter)]
only_m=m_id[~np.in1d(m_id,inter)]

u_ms=stms[np.in1d(subfindID,union)]
print u_ms[:10]
in_ms=stms[np.in1d(subfindID,inter)]

e_ms=stms[~np.in1d(subfindID,union)]

print u_ms.shape,in_ms.shape,e_ms.shape,tri_id.shape,only_k.shape,only_m.shape




## number counts
plt.xlabel('log10(stellar mass)')
plt.ylabel('number counts')
plt.title('Central non-disk galaxy')
plt.hist(np.log10(e_ms),bins=20)
plt.savefig(dirc+'compare/hist_099_eC.png')
#plt.show()

