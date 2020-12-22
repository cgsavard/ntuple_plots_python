import uproot
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep

def res_plot(arrays, x_feat, res_feat):
    # cut where there are matches
    tp_cut = arrays["tp_nmatch"].flatten()>0
    
    # select feature being used and create bins
    if x_feat == "pt":
        tp_feat = arrays["tp_pt"].flatten()[tp_cut]
        feat_bins = np.logspace(.3,2,20)
        feat_dig = np.digitize(tp_feat,feat_bins)
    elif x_feat == "eta":
        tp_feat = arrays["tp_eta"].flatten()[tp_cut]
        feat_bins = np.linspace(-2.4,2.4,30)
        feat_dig = np.digitize(tp_feat,feat_bins)
    elif x_feat == "phi":
        tp_feat = arrays["tp_phi"].flatten()[tp_cut]
        feat_bins = np.linspace(-3.2,3.2,30)
        feat_dig = np.digitize(tp_feat,feat_bins)
    elif x_feat == "d0":
        tp_feat = abs(arrays["tp_d0"].flatten())[tp_cut]
        feat_bins = np.linspace(0,1,10)
        feat_dig = np.digitize(tp_feat,feat_bins)
    elif x_feat == "z0":
        tp_feat = arrays["tp_z0"].flatten()[tp_cut]
        feat_bins = np.linspace(-15,15,30)
        feat_dig = np.digitize(tp_feat,feat_bins)
    else:
        print("not valid feature")
        return None
      
    # loop through each bin and calculate resolution and error
    res = np.zeros(len(feat_bins)-1)
    res_err = np.zeros(len(feat_bins)-1)
    for i in range(1,len(feat_bins)):
        if res_feat == "pt weighted":
            trk_res_temp = (abs(arrays["tp_pt"].flatten()[tp_cut]-arrays["matchtrk_pt"].flatten()[tp_cut])/(arrays["tp_pt"].flatten()[tp_cut]))[feat_dig==i]
        else:
            trk_res_temp = abs(arrays["tp_"+res_feat].flatten()[tp_cut]-arrays["matchtrk_"+res_feat].flatten()[tp_cut])[feat_dig==i]
        
        res[i-1] = np.mean(trk_res_temp)
        res_err[i-1] = np.std(trk_res_temp)/np.sqrt(len(trk_res_temp))
    
    # correct bin to be centered and error is bin width
    bins = [(feat_bins[i]+feat_bins[i+1])/2 for i in range(len(feat_bins)-1)]
    bins_err = [(feat_bins[i+1]-feat_bins[i])/2 for i in range(len(feat_bins)-1)]
    
    return bins, bins_err, res, res_err


# take in .root into arrays
arrays_disp = (uproot.open("../data/TTbar_pu200_D49_extended.root")["L1TrackNtuple/eventTree"].arrays("*", namedecode="utf-8"))
arrays_pro = (uproot.open("../data/TTbar_pu200_D49_prompt.root")["L1TrackNtuple/eventTree"].arrays("*", namedecode="utf-8"))

# plot 2D resolution
bins, bins_err, res_disp, res_err_disp = res_plot(arrays_disp, "eta","eta")
bins, bins_err, res_pro, res_err_pro = res_plot(arrays_pro, "eta","eta")
plt.style.use(hep.style.CMS)
plt.errorbar(bins,res_disp,xerr=bins_err,yerr=res_err_disp,fmt='o',label='ttbar displaced')
plt.errorbar(bins,res_pro,xerr=bins_err,yerr=res_err_pro,fmt='o',label='ttbar prompt')
plt.legend(loc="best")
plt.xlabel(r"tracking particle $\eta$")
plt.ylabel(r"|$\eta$ resolution|")
plt.show()

bins, bins_err, res_disp, res_err_disp = res_plot(arrays_disp, "eta","pt")
bins, bins_err, res_pro, res_err_pro = res_plot(arrays_pro, "eta","pt")
plt.style.use(hep.style.CMS)
plt.errorbar(bins,res_disp,xerr=bins_err,yerr=res_err_disp,fmt='o',label='ttbar displaced')
plt.errorbar(bins,res_pro,xerr=bins_err,yerr=res_err_pro,fmt='o',label='ttbar prompt')
plt.legend(loc="best")
plt.xlabel(r"tracking particle $\eta$")
plt.ylabel(r"|$p_T$ resolution|")
plt.show()

bins, bins_err, res_disp, res_err_disp = res_plot(arrays_disp, "eta","pt weighted")
bins, bins_err, res_pro, res_err_pro = res_plot(arrays_pro, "eta","pt weighted")
plt.style.use(hep.style.CMS)
plt.errorbar(bins,res_disp,xerr=bins_err,yerr=res_err_disp,fmt='o',label='ttbar displaced')
plt.errorbar(bins,res_pro,xerr=bins_err,yerr=res_err_pro,fmt='o',label='ttbar prompt')
plt.legend(loc="best")
plt.xlabel(r"tracking particle $\eta$")
plt.ylabel(r"|$p_T$ resolution| / tracking particle $p_T$")
plt.show()

bins, bins_err, res_disp, res_err_disp = res_plot(arrays_disp, "eta","phi")
bins, bins_err, res_pro, res_err_pro = res_plot(arrays_pro, "eta","phi")
plt.style.use(hep.style.CMS)
plt.errorbar(bins,res_disp,xerr=bins_err,yerr=res_err_disp,fmt='o',label='ttbar displaced')
plt.errorbar(bins,res_pro,xerr=bins_err,yerr=res_err_pro,fmt='o',label='ttbar prompt')
plt.legend(loc="best")
plt.xlabel(r"tracking particle $\eta$")
plt.ylabel(r"|$\phi$ resolution|")
plt.show()
