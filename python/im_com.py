import matplotlib.pyplot as plt
import matplotlib.image as im

#im1=im.imread('/Users/jwhsueh/Documents/SHARP_jw/data/glamer/glamer_tri_fold_pd.png')
#im2=im.imread('/Users/jwhsueh/Documents/SHARP_jw/data/glamer/glamer_elp_fold_pd.png')
#im3=im.imread('/Users/jwhsueh/Documents/SHARP_jw/data/glamer/glamer_edg_fold_pd.png')
#im4=im.imread('/Users/jwhsueh/Documents/SHARP_jw/data/glamer/glamer_fa_fold_pd.png')

#im1=im.imread('/Users/jwhsueh/Documents/SHARP_jw/data/glamer/glamer_sie_cusp_pd.png')
#im2=im.imread('/Users/jwhsueh/Documents/SHARP_jw/data/glamer/glamer_sie_fold_pd.png')

im1=im.imread('../data/VLBA_b0712_1.ps')

#im1=im.imread('/Users/jwhsueh/Documents/SHARP_jw/data/sub_gravlens/snap99_ast/ast_anomaly_cusp.png')
#im2=im.imread('/Users/jwhsueh/Documents/SHARP_jw/data/sub_gravlens/snap99_ast/ast_anomaly_fold.png')
'''
fig, ((ax1,ax2)) = plt.subplots(1,2)

ax1.imshow(im1)
ax1.axis('off')
ax2.imshow(im2)
ax2.axis('off')
#ax1.imshow(im3)
#ax1.axis('off')
#ax2.imshow(im4)
#ax2.axis('off')

#fig.text(0.04,0.8,'i=7',color='w')
plt.subplots_adjust(wspace=0,hspace=0)
'''
plt.imshow(im1)
plt.show()
#plt.savefig('/Users/jwhsueh/Documents/SHARP_jw/data/glamer/glamer_sie_pd.png',bbox_inches='tight',dpi=200)
