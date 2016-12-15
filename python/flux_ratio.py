import numpy as np
import matplotlib.pyplot as plt

mr_15ghz=np.array([15.6,11.6,5.9,1.1])
mr_15ghz_sig = mr_15ghz*0.2
mr_15ghz_sig = np.array([0.78,0.58,0.3,0.06])
mr_15ghz_mod = mr_15ghz_sig/mr_15ghz

vlba_165 = np.array([10.1,6.7,3.3,0.37])
vlba_165_sig = vlba_165*0.2
vlba_165_sig = np.array([0.5,0.3,0.2,0.01])
vlba_165_mod = vlba_165_sig/vlba_165

mr_5ghz_fr = np.array([0.843,0.418,0.082])
mr_5ghz_fr_sig = np.array([0.061,0.037,0.035])

mr_15ghz_fr = np.array([mr_15ghz[1]/mr_15ghz[0],mr_15ghz[2]/mr_15ghz[0],mr_15ghz[3]/mr_15ghz[0]])
mr_15ghz_fr_sig = mr_15ghz_fr*0.4
#mr_15ghz_fr_sig.fill(0.283)
mr_15ghz_fr_sig = np.array([mr_15ghz_fr[0]*(mr_15ghz_mod[0]+mr_15ghz_mod[1]),mr_15ghz_fr[1]*(mr_15ghz_mod[0]+mr_15ghz_mod[2]),mr_15ghz_fr[2]*(mr_15ghz_mod[0]+mr_15ghz_mod[3])])
print mr_15ghz_fr
print mr_15ghz_fr_sig

vlba_165_fr = np.array([vlba_165[1]/vlba_165[0],vlba_165[2]/vlba_165[0],vlba_165[3]/vlba_165[0]])
vlba_165_fr_sig = vlba_165_fr*0.4
#vlba_165_fr_sig.fill(0.283)
vlba_165_fr_sig = np.array([vlba_165_fr[0]*(vlba_165_mod[0]+vlba_165_mod[1]),vlba_165_fr[1]*(vlba_165_mod[0]+vlba_165_mod[2]),vlba_165_fr[2]*(vlba_165_mod[0]+vlba_165_mod[3])])
print vlba_165_fr
print vlba_165_fr_sig

fig,(fig1,fig2,fig3)=plt.subplots(3,sharex=True)
## B/A
fig1.errorbar(15,mr_15ghz_fr[0],yerr=(mr_15ghz_fr_sig[0]),color='b')
fig1.errorbar(5,mr_5ghz_fr[0],yerr=(mr_5ghz_fr_sig[0]),color='r')
fig1.errorbar(1.65,vlba_165_fr[0],yerr=(vlba_165_fr_sig[0]),color='k')
fig1.scatter(15,mr_15ghz_fr[0],color='b',marker='^')
fig1.scatter(5,mr_5ghz_fr[0],color='r',marker='^')
fig1.scatter(1.65,vlba_165_fr[0],color='k',marker='^')
fig1.set_ylabel('(B/A)')

## C/A
fig2.errorbar(15,mr_15ghz_fr[1],yerr=(mr_15ghz_fr_sig[1]),color='b')
fig2.errorbar(5,mr_5ghz_fr[1],yerr=(mr_5ghz_fr_sig[1]),color='r')
fig2.errorbar(1.65,vlba_165_fr[1],yerr=(vlba_165_fr_sig[1]),color='k')
fig2.scatter(15,mr_15ghz_fr[1],color='b',marker='^')
fig2.scatter(5,mr_5ghz_fr[1],color='r',marker='^')
fig2.scatter(1.65,vlba_165_fr[1],color='k',marker='^')
fig2.set_ylabel('(C/A)')

## D/A
fig3.errorbar(15,mr_15ghz_fr[2],yerr=(mr_15ghz_fr_sig[2]),color='b')
fig3.errorbar(5,mr_5ghz_fr[2],yerr=(mr_5ghz_fr_sig[2]),color='r')
fig3.errorbar(1.65,vlba_165_fr[2],yerr=(vlba_165_fr_sig[2]),color='k')
fig3.scatter(15,mr_15ghz_fr[2],color='b',marker='^')
fig3.scatter(5,mr_5ghz_fr[2],color='r',marker='^')
fig3.scatter(1.65,vlba_165_fr[2],color='k',marker='^')
fig3.set_ylabel('(D/A)')

plt.xlabel('Frequency (GHz)')
#fig.set_ylabel('Flux-ratio density')

#fig.subplots_adjust(hspace=0)
#plt.show()

plt.savefig('../data/B0712_flux_ratio.png')