import numpy as np
import groupcat

dirc='../../data/illustris_1/'

snapshotPath = '/Volumes/narsil_1/jwhsueh/illustris_1'
snapNum = '099'

class cosmopara:
	h = 0.704
	OM = 0.27

## drawing number

#e_num=50
#tri_num=22
#dou_num=88
#only_m_num=44
#only_k_num=66

e_num=65
tri_num=0
dou_num=0
only_m_num=0
only_k_num=0


## all
table=np.loadtxt(dirc+'AllGalaxy_099_x.dat',skiprows=1)
table2=np.loadtxt(dirc+'AllGalaxy_099_y.dat',skiprows=1)
table3=np.loadtxt(dirc+'AllGalaxy_099_z.dat',skiprows=1)

table4=np.loadtxt(dirc+'Galaxy_99_cen.dat')

### central galaxy

cenID=table4[:,0]
subDM=table4[:,5]

# DM mass cut
mask=subDM<5e13
cenID=cenID[mask]

## ------------------- ##

## projection1
subfindID=table[:,0]
df=table[:,1]
bf=table[:,2]
Re=table[:,3]
#subflag=table[:,5]
ms=table[:,6]
ser=table[:,8]
mor=table[:,10]

##--- central & DM mass cut
mask=np.in1d(subfindID,cenID)
subfindID,bf,ser,mor,df,Re,ms=subfindID[mask],bf[mask],ser[mask],mor[mask],df[mask],Re[mask],ms[mask]
Rea,msa = Re,ms

all_list = '../../data/illustris_1/snap99_p1_all.txt'
np.savetxt(all_list,np.c_[msa,Rea],newline='\n')

## Re flag
#mask=Re>0.3
#subfindID,bf,ser,mor,df=subfindID[mask],bf[mask],ser[mask],mor[mask],df[mask]

## set up subgroup

bf_mask=bf<0.6
df_mask=df>0.4
mor_mask=(ser<2.0)&(mor==0.)
tri_mask=(bf<0.6)&(df>0.4)&(ser<2.0)&(mor==0.)

k_id=subfindID[bf_mask]
k_ms=ms[bf_mask]
m_id=subfindID[mor_mask]
m_ms=ms[mor_mask] 
m_re=Re[mor_mask]

tri_id=subfindID[tri_mask]
tri_ms=ms[tri_mask]
tri_re=Re[tri_mask]
tri_bf=bf[tri_mask]
tri_df=df[tri_mask]
tri_dim = tri_ms*(1-2.*tri_bf)

d_id = subfindID[df_mask]
d_ms=ms[df_mask]
d_re=Re[df_mask]

mor_file = '../../data/illustris_1/snap99_p1_mor.txt'
np.savetxt(mor_file,np.c_[m_ms,m_re])
k_file = '../../data/illustris_1/snap99_p1_k.txt'
np.savetxt(k_file,np.c_[d_ms,d_re])
tri_file = '../../data/illustris_1/snap99_p1_tri.txt'
np.savetxt(tri_file,np.c_[tri_id,tri_ms,tri_re,tri_dim,tri_df])

u_id=np.union1d(k_id,m_id)
in_id=np.intersect1d(k_id,m_id)
dou_id=in_id[~np.in1d(in_id,tri_id)]

only_k=k_id[~np.in1d(k_id,in_id)]
only_m=m_id[~np.in1d(m_id,in_id)]

e_mask = ~np.in1d(subfindID,u_id)
e_id=subfindID[e_mask]
e_ms=ms[e_mask]
e_re = Re[e_mask]

e_file = '../../data/illustris_1/snap99_p1_e.txt'
np.savetxt(e_file,np.c_[e_ms,e_re])

#print e_id.shape,tri_id.shape,dou_id.shape,only_k.shape,only_m.shape

# projection1 mark
e_idx,tri_idx,dou_idx,only_mx,only_kx=e_id+0.1,tri_id+0.1,dou_id+0.1,only_m+0.1,only_k+0.1

## -------------  ##

## projection2
subfindID=table2[:,0]
df=table2[:,1]
bf=table2[:,2]
Re=table2[:,3]
#subflag=table[:,5]
ms=table[:,6]
ser=table2[:,8]
mor=table2[:,10]

##--- central & DM mass cut
mask=np.in1d(subfindID,cenID)
subfindID,bf,ser,mor,df,Re,ms=subfindID[mask],bf[mask],ser[mask],mor[mask],df[mask],Re[mask],ms[mask]
Rea,msa = Re,ms

all_list = '../../data/illustris_1/snap99_p2_all.txt'
np.savetxt(all_list,np.c_[msa,Rea],newline='\n')

## Re flag
#mask=Re>0.0
#subfindID,bf,ser,mor,df=subfindID[mask],bf[mask],ser[mask],mor[mask],df[mask]

## set up subgroup

bf_mask=bf<0.6
df_mask=df>0.4
mor_mask=(ser<2.0)&(mor==0.)
tri_mask=(bf<0.6)&(df>0.4)&(ser<2.0)&(mor==0.)

k_id=subfindID[bf_mask]
k_ms=ms[bf_mask]
m_id=subfindID[mor_mask]
m_ms=ms[mor_mask] 
m_re=Re[mor_mask]

tri_id=subfindID[tri_mask]
tri_ms=ms[tri_mask]
tri_re=Re[tri_mask]
tri_bf=bf[tri_mask]
tri_df=df[tri_mask]
tri_dim = tri_ms*(1-2.*tri_bf)

d_id = subfindID[df_mask]
d_ms=ms[df_mask]
d_re=Re[df_mask]

mor_file = '../../data/illustris_1/snap99_p2_mor.txt'
np.savetxt(mor_file,np.c_[m_ms,m_re])
k_file = '../../data/illustris_1/snap99_p2_k.txt'
np.savetxt(k_file,np.c_[d_ms,d_re])
tri_file = '../../data/illustris_1/snap99_p2_tri.txt'
np.savetxt(tri_file,np.c_[tri_id,tri_ms,tri_re,tri_dim,tri_df])

u_id=np.union1d(k_id,m_id)
in_id=np.intersect1d(k_id,m_id)
dou_id=in_id[~np.in1d(in_id,tri_id)]

only_k=k_id[~np.in1d(k_id,in_id)]
only_m=m_id[~np.in1d(m_id,in_id)]

e_mask = ~np.in1d(subfindID,u_id)
e_id=subfindID[e_mask]
e_ms=ms[e_mask]
e_re = Re[e_mask]

e_file = '../../data/illustris_1/snap99_p2_e.txt'
np.savetxt(e_file,np.c_[e_ms,e_re])

#print e_id.shape,tri_id.shape,dou_id.shape,only_k.shape,only_m.shape

# projection2 mark
e_idy,tri_idy,dou_idy,only_my,only_ky=e_id+0.2,tri_id+0.2,dou_id+0.2,only_m+0.2,only_k+0.2

## --------------- ##

## projection3
subfindID=table3[:,0]
df=table3[:,1]
bf=table3[:,2]
Re=table3[:,3]
#subflag=table[:,5]
ms=table[:,6]
ser=table3[:,8]
mor=table3[:,10]

##--- central & DM mass cut
mask=np.in1d(subfindID,cenID)
subfindID,bf,ser,mor,df,Re,ms=subfindID[mask],bf[mask],ser[mask],mor[mask],df[mask],Re[mask],ms[mask]
Rea,msa = Re,ms

all_list = '../../data/illustris_1/snap99_p3_all.txt'
np.savetxt(all_list,np.c_[msa,Rea],newline='\n')

## Re flag
#mask=Re>0.0
#subfindID,bf,ser,mor,df=subfindID[mask],bf[mask],ser[mask],mor[mask],df[mask]

## set up subgroup

bf_mask=bf<0.6
df_mask=df>0.4
mor_mask=(ser<2.0)&(mor==0.)
tri_mask=(bf<0.6)&(df>0.4)&(ser<2.0)&(mor==0.)

k_id=subfindID[bf_mask]
k_ms=ms[bf_mask]
m_id=subfindID[mor_mask]
m_ms=ms[mor_mask] 
m_re=Re[mor_mask]

tri_id=subfindID[tri_mask]
tri_ms=ms[tri_mask]
tri_re=Re[tri_mask]
tri_bf=bf[tri_mask]
tri_df=df[tri_mask]
tri_dim = tri_ms*(1-2.*tri_bf)

d_id = subfindID[df_mask]
d_ms=ms[df_mask]
d_re=Re[df_mask]

mor_file = '../../data/illustris_1/snap99_p3_mor.txt'
np.savetxt(mor_file,np.c_[m_ms,m_re])
k_file = '../../data/illustris_1/snap99_p3_k.txt'
np.savetxt(k_file,np.c_[d_ms,d_re])
tri_file = '../../data/illustris_1/snap99_p3_tri.txt'
np.savetxt(tri_file,np.c_[tri_id,tri_ms,tri_re,tri_dim,tri_df])

u_id=np.union1d(k_id,m_id)
in_id=np.intersect1d(k_id,m_id)
dou_id=in_id[~np.in1d(in_id,tri_id)]

only_k=k_id[~np.in1d(k_id,in_id)]
only_m=m_id[~np.in1d(m_id,in_id)]

e_mask = ~np.in1d(subfindID,u_id)
e_id=subfindID[e_mask]
e_ms=ms[e_mask]
e_re = Re[e_mask]

e_file = '../../data/illustris_1/snap99_p3_e.txt'
np.savetxt(e_file,np.c_[e_ms,e_re])

#print e_id.shape,tri_id.shape,dou_id.shape,only_k.shape,only_m.shape

# projection3 mark
e_idz,tri_idz,dou_idz,only_mz,only_kz=e_id+0.3,tri_id+0.3,dou_id+0.3,only_m+0.3,only_k+0.3


## ----------- ##



e_id=np.append(e_idx,e_idy)
e_id=np.append(e_id,e_idz)

tri_id=np.append(tri_idx,tri_idy)
tri_id=np.append(tri_id,tri_idz)

dou_id=np.append(dou_idx,dou_idy)
dou_id=np.append(dou_id,dou_idz)

only_m=np.append(only_mx,only_my)
only_m=np.append(only_m,only_mz)

only_k=np.append(only_kx,only_ky)
only_k=np.append(only_k,only_kz)

print tri_id.size,tri_id

## start drawing

e_draw=np.sort(np.random.choice(e_id,e_num,replace=False))
#print e_draw
tri_draw=np.sort(np.random.choice(tri_id,tri_num,replace=False))
dou_draw=np.sort(np.random.choice(dou_id,dou_num,replace=False))
only_m_draw=np.sort(np.random.choice(only_m,only_m_num,replace=False))
only_k_draw=np.sort(np.random.choice(only_k,only_k_num,replace=False))

## --------------- ##

## create catalog

# from drawing
e_int=e_draw.astype(int)
tri_int=tri_draw.astype(int)
dou_int=dou_draw.astype(int)
only_m_int=only_m_draw.astype(int)
only_k_int=only_k_draw.astype(int)

catalogID=np.union1d(e_int,tri_int)
catalogID=np.union1d(catalogID,dou_int)
catalogID=np.union1d(catalogID,only_m_int)
catalogID=np.union1d(catalogID,only_k_int)

cat_full=np.empty((len(catalogID),4))
cat_full[:,0]=catalogID
cat_full[:,1:].fill(1) # fill with False flag

projID=np.append(e_draw,tri_draw)
projID=np.append(projID,dou_draw)
projID=np.append(projID,only_m_draw)
projID=np.append(projID,only_k_draw)

projID=np.sort(projID)
print projID[:20]

for element in projID:
	galaxyID=np.floor(element)

	proj=np.round((element-galaxyID)*10.)

	idx=list(catalogID).index(galaxyID)

	cat_full[idx,proj]=0

cat_full=cat_full.astype(int)
print cat_full[:10,:]

#np.savetxt(dirc+'raytrace_catalog2_elp.dat',cat_full,fmt='%d',delimiter='\t')

## ---------- ##

e_projID=e_draw
d_projID=np.append(tri_draw,dou_draw)
d_projID=np.append(d_projID,only_m_draw)
d_projID=np.append(d_projID,only_k_draw)

cat_mark=np.empty((len(catalogID),4))
cat_mark[:,0]=catalogID
cat_mark[:,1:].fill(1) # fill with False flag

for element in e_projID:
	galaxyID=np.floor(element)
	
	proj=np.round((element-galaxyID)*10.)

	idx=list(catalogID).index(galaxyID)
	cat_mark[idx,proj]=3

for element in d_projID:
	galaxyID=np.floor(element)
	proj=np.round((element-galaxyID)*10.)
	idx=list(catalogID).index(galaxyID)
	cat_mark[idx,proj]=2

#np.savetxt(dirc+'raytrace_catalog2_tri.dat',cat_mark,fmt='%d',delimiter='\t')
