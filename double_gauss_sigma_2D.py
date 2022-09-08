import uproot3
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
from sklearn.mixture import GaussianMixture

plt.style.use(hep.style.CMS)

# take in .root into arrays
data = (uproot3.open("../2022_code/d0_bias_fix/DispMu_PU200_D49_extended_before.root")["L1TrackNtuple/eventTree"].arrays("*",namedecode="utf-8"))

def sig_plot(arrays, x_feat, res_feat, ran=1):
    # cut where there are matches
    tp_cut = arrays["tp_nmatch"].flatten()>0
    #tp_cut = np.logical_and(arrays["tp_nmatch"].flatten()>0,abs(arrays['tp_eta'].flatten())>1.8)
    
    # select feature being used and create bins
    if x_feat == "pt":
        tp_feat = arrays["tp_pt"].flatten()
        feat_bins = [2,3,5,10,50]
        feat_dig = np.digitize(tp_feat,feat_bins)
    elif x_feat == "eta":
        tp_feat = arrays["tp_eta"].flatten()
        feat_bins = [0,.7,1,1.6,2.4]
        feat_dig = np.digitize(abs(tp_feat),feat_bins)
    elif x_feat == "phi":
        tp_feat = arrays["tp_phi"].flatten()
        feat_bins = [0,1,2,3.4]
        feat_dig = np.digitize(abs(tp_feat),feat_bins)
    elif x_feat == "d0":
        tp_feat = arrays["tp_d0"].flatten()
        feat_bins = [-10,-4,-2,-1,-.5,0,0.5,1,2,4,10]
        feat_dig = np.digitize(tp_feat,feat_bins)
    elif x_feat == "z0":
        tp_feat = arrays["tp_z0"].flatten()
        feat_bins = [0,5,10,15,20,30]
        feat_dig = np.digitize(abs(tp_feat),feat_bins)
    else:
        print("not valid feature")
        return None
      
    # loop through each bin and calculate resolution and error
    sig = np.zeros(len(feat_bins)-1)
    sig_err = np.zeros(len(feat_bins)-1)
    for i in range(1,len(feat_bins)):
        
        res = (arrays["tp_"+res_feat].flatten()-arrays["matchtrk_"+res_feat].flatten())[np.logical_and(tp_cut,feat_dig==i)].reshape(-1, 1)
        res = res[np.isfinite(res)].reshape(-1, 1)
        if len(res)==0:
            sig[i-1] = -1
            continue
        gm = GaussianMixture(n_components=2, random_state=0).fit(res)
        
        textstr = '\n'.join((r'$%.2f\leq$tp_%s$<%.2f$' % (feat_bins[i-1], x_feat, feat_bins[i]),
                         r'%d entries' % (len(res),),
                         r'$\mu_1=%.2f$, $\mu_2=%.2f$' % (gm.means_[0], gm.means_[1]),
                         r'$\sigma_1=%.2f$, $\sigma_2=%.2f$' % (np.sqrt(gm.covariances_[0]), np.sqrt(gm.covariances_[1])),
                         r'$w_1=%.3f$, $w_2=%.3f$' % (gm.weights_[0], gm.weights_[1])))
        #print(textstr)
        sig[i-1] = np.sqrt(gm.weights_[0]*gm.covariances_[0]+gm.weights_[1]*gm.covariances_[1])
        sig_err[i-1] = 0 # for now, don't have errors
        
        # plot original and double gauss to check fit
        t = np.linspace(-ran, ran, 100).reshape(100,1)
        logprob = gm.score_samples(t)
        pdf = np.exp(logprob)
        plt.hist(res,bins=100,range=(-ran,ran),histtype='step',density=True)
        plt.plot(t,pdf,'-k')
        plt.xlabel(x_feat+ ' resolution')
        plt.ylabel('normalized counts')
        plt.text(0.75,.85, textstr,horizontalalignment='center',verticalalignment='center',transform=plt.gca().transAxes)
        plt.show()
    
    # correct bin to be centered and error is bin width
    bins = [(feat_bins[i]+feat_bins[i+1])/2 for i in range(len(feat_bins)-1)]
    bins_err = [(feat_bins[i+1]-feat_bins[i])/2 for i in range(len(feat_bins)-1)]
    
    return bins, bins_err, sig, sig_err


bins, bins_err, sig, sig_err = sig_plot(data, "d0","d0",1)
plt.errorbar(bins,sig,xerr=bins_err,yerr=sig_err,fmt='o')
plt.xlabel(r"tracking particle $d_0$ (cm)")
plt.ylabel(r"$\sigma_{d_0 res}$ (cm)")
plt.ylim(bottom=0)
plt.show()
