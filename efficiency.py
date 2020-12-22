import uproot
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep

def eff_plot(arrays, feat):
    # tp quality cuts  
    tp_cuts = np.logical_and(np.logical_and(np.logical_and(arrays["tp_pt"].flatten()>=2,
                                                           arrays["tp_pt"].flatten()<100),
                                                           abs(arrays["tp_eta"].flatten())<2.4),
                                                           arrays["tp_eventid"].flatten()==0) # only look at hard interaction
    
    # select feature being used and create bins
    if feat == "pt":
        tp_feat = arrays["tp_pt"].flatten()[tp_cuts]
        feat_bins = np.logspace(.3,2,20)
        feat_dig = np.digitize(tp_feat,feat_bins)
    elif feat == "eta":
        tp_feat = arrays["tp_eta"].flatten()[tp_cuts]
        feat_bins = np.linspace(-2.4,2.4,30)
        feat_dig = np.digitize(tp_feat,feat_bins)
    elif feat == "phi":
        tp_feat = arrays["tp_phi"].flatten()[tp_cuts]
        feat_bins = np.linspace(-3.2,3.2,30)
        feat_dig = np.digitize(tp_feat,feat_bins)
    elif feat == "d0":
        tp_feat = abs(arrays["tp_d0"].flatten()[tp_cuts])
        feat_bins = np.linspace(0,1,10)
        feat_dig = np.digitize(tp_feat,feat_bins)
    elif feat == "z0":
        tp_feat = arrays["tp_z0"].flatten()[tp_cuts]
        feat_bins = np.linspace(-15,15,30)
        feat_dig = np.digitize(tp_feat,feat_bins)
    elif feat == "dxy":
        tp_feat = arrays["tp_dxy"].flatten()[tp_cuts]
        feat_bins = np.linspace(0,50,10)
        feat_dig = np.digitize(tp_feat,feat_bins)
    else:
        print("not valid feature")
        return None
    
    # find the matching of tps
    tp_nmatch = arrays["tp_nmatch"].flatten()[tp_cuts]
    
    # loop through each bin and calculate efficiency and error
    eff = np.zeros(len(feat_bins)-1)
    eff_err = np.zeros(len(feat_bins)-1)
    for i in range(1,len(feat_bins)):
        tp_nmatch_temp = tp_nmatch[feat_dig==i]
        if len(tp_nmatch_temp)==0:
            continue
        
        eff[i-1] = np.sum(tp_nmatch_temp>0)/len(tp_nmatch_temp)
        eff_err[i-1] = np.sqrt(eff[i-1]*(1-eff[i-1])/len(tp_nmatch_temp))
    
    # correct bin to be centered and error is bin width
    bins = [(feat_bins[i]+feat_bins[i+1])/2 for i in range(len(feat_bins)-1)]
    bins_err = [(feat_bins[i+1]-feat_bins[i])/2 for i in range(len(feat_bins)-1)]
    
    return bins, bins_err, eff, eff_err


# take in .root into arrays
arrays_disp = (uproot.open("../data/TTbar_pu200_D49_extended.root")["L1TrackNtuple/eventTree"].arrays("*", namedecode="utf-8"))
arrays_pro = (uproot.open("../data/TTbar_pu200_D49_prompt.root")["L1TrackNtuple/eventTree"].arrays("*", namedecode="utf-8"))

# efficiency plots
bins, bins_err, eff_disp, eff_err_disp = eff_plot(arrays_disp, "pt")
bins, bins_err, eff_pro, eff_err_pro = eff_plot(arrays_pro, "pt")
plt.style.use(hep.style.CMS)
plt.errorbar(bins,eff_disp,xerr=bins_err,yerr=eff_err_disp,fmt='o',label='ttbar displaced')
plt.errorbar(bins,eff_pro,xerr=bins_err,yerr=eff_err_pro,fmt='o',label='ttbar prompt')
plt.legend(loc="best")
plt.xlabel(r"tracking particle $p_T$ (GeV/c)")
plt.ylabel("efficiency")
plt.ylim(0,1.1)
plt.show()

bins, bins_err, eff_disp, eff_err_disp = eff_plot(arrays_disp, "eta")
bins, bins_err, eff_pro, eff_err_pro = eff_plot(arrays_pro, "eta")
plt.style.use(hep.style.CMS)
plt.errorbar(bins,eff_disp,xerr=bins_err,yerr=eff_err_disp,fmt='o',label='ttbar displaced')
plt.errorbar(bins,eff_pro,xerr=bins_err,yerr=eff_err_pro,fmt='o',label='ttbar prompt')
plt.legend(loc="best")
plt.xlabel(r"tracking particle $\eta$")
plt.ylabel("efficiency")
plt.ylim(0,1.1)
plt.show()

bins, bins_err, eff_disp, eff_err_disp = eff_plot(arrays_disp, "phi")
bins, bins_err, eff_pro, eff_err_pro = eff_plot(arrays_pro, "phi")
plt.style.use(hep.style.CMS)
plt.errorbar(bins,eff_disp,xerr=bins_err,yerr=eff_err_disp,fmt='o',label='ttbar displaced')
plt.errorbar(bins,eff_pro,xerr=bins_err,yerr=eff_err_pro,fmt='o',label='ttbar prompt')
plt.legend(loc="best")
plt.xlabel(r"tracking particle $\phi$")
plt.ylabel("efficiency")
plt.ylim(0,1.1)
plt.show()

bins, bins_err, eff_disp, eff_err_disp = eff_plot(arrays_disp, "z0")
bins, bins_err, eff_pro, eff_err_pro = eff_plot(arrays_pro, "z0")
plt.style.use(hep.style.CMS)
plt.errorbar(bins,eff_disp,xerr=bins_err,yerr=eff_err_disp,fmt='o',label='ttbar displaced')
plt.errorbar(bins,eff_pro,xerr=bins_err,yerr=eff_err_pro,fmt='o',label='ttbar prompt')
plt.legend(loc="best")
plt.xlabel(r"tracking particle $z_0$ (cm)")
plt.ylabel("efficiency")
plt.ylim(0,1.1)
plt.show()

bins, bins_err, eff_disp, eff_err_disp = eff_plot(arrays_disp, "d0")
bins, bins_err, eff_pro, eff_err_pro = eff_plot(arrays_pro, "d0")
plt.style.use(hep.style.CMS)
plt.errorbar(bins,eff_disp,xerr=bins_err,yerr=eff_err_disp,fmt='o',label='ttbar displaced')
plt.errorbar(bins,eff_pro,xerr=bins_err,yerr=eff_err_pro,fmt='o',label='ttbar prompt')
plt.legend(loc="best")
plt.xlabel(r"tracking particle |$d_0$| (cm)")
plt.ylabel("efficiency")
plt.ylim(0,1.1)
plt.show()

bins, bins_err, eff_disp, eff_err_disp = eff_plot(arrays_disp, "dxy")
bins, bins_err, eff_pro, eff_err_pro = eff_plot(arrays_pro, "dxy")
plt.style.use(hep.style.CMS)
plt.errorbar(bins,eff_disp,xerr=bins_err,yerr=eff_err_disp,fmt='o',label='ttbar displaced')
plt.errorbar(bins,eff_pro,xerr=bins_err,yerr=eff_err_pro,fmt='o',label='ttbar prompt')
plt.legend(loc="best")
plt.xlabel(r"tracking particle |$d_{xy}$| (cm)")
plt.ylabel("efficiency")
plt.ylim(0,1.1)
plt.show()
