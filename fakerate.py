import uproot
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep

def faker_plot(arrays, feat):
    # trk quality cuts
    trk_cuts = np.logical_and(np.logical_and(arrays["trk_pt"].flatten()>=2,
                                             arrays["trk_nstub"].flatten()>=4),
                                             abs(arrays["trk_eta"].flatten())<2.4)
    
    # select feature being used and create bins
    if feat == "pt":
        trk_feat = arrays["trk_pt"].flatten()[trk_cuts]
        feat_bins = np.logspace(.3,2,20)
        feat_dig = np.digitize(trk_feat,feat_bins)
    elif feat == "eta":
        trk_feat = arrays["trk_eta"].flatten()[trk_cuts]
        feat_bins = np.linspace(-2.4,2.4,30)
        feat_dig = np.digitize(trk_feat,feat_bins)
    elif feat == "phi":
        trk_feat = arrays["trk_phi"].flatten()[trk_cuts]
        feat_bins = np.linspace(-3.2,3.2,30)
        feat_dig = np.digitize(trk_feat,feat_bins)
    elif feat == "z0":
        trk_feat = arrays["trk_z0"].flatten()[trk_cuts]
        feat_bins = np.linspace(-15,15,30)
        feat_dig = np.digitize(trk_feat,feat_bins)
    elif feat == "d0":
        trk_feat = abs(arrays["trk_d0"].flatten()[trk_cuts])
        feat_bins = np.linspace(0,1,30)
        feat_dig = np.digitize(trk_feat,feat_bins)
    elif feat == "chi2":
        trk_feat = arrays["trk_chi2"].flatten()[trk_cuts]
        feat_bins = np.linspace(0,60,30)
        feat_dig = np.digitize(trk_feat,feat_bins)
    else:
        print("not valid feature")
        return None
    
    # find the matching of tps
    trk_fake = arrays["trk_fake"].flatten()[trk_cuts]
    
    # loop through each bin and calculate efficiency and error
    faker = np.zeros(len(feat_bins)-1)
    faker_err = np.zeros(len(feat_bins)-1)
    for i in range(1,len(feat_bins)):
        trk_fake_temp = trk_fake[feat_dig==i]
        if len(trk_fake_temp)==0:
            continue
        
        faker[i-1] = np.sum(trk_fake_temp==0)/len(trk_fake_temp)
        faker_err[i-1] = np.sqrt(faker[i-1]*(1-faker[i-1])/len(trk_fake_temp))
    
    # correct bin to be centered and error is bin width
    bins = [(feat_bins[i]+feat_bins[i+1])/2 for i in range(len(feat_bins)-1)]
    bins_err = [(feat_bins[i+1]-feat_bins[i])/2 for i in range(len(feat_bins)-1)]
    
    return bins, bins_err, faker, faker_err


# take in .root into arrays
arrays_disp = (uproot.open("../data/TTbar_pu200_D49_extended.root")["L1TrackNtuple/eventTree"].arrays("*", namedecode="utf-8"))
arrays_pro = (uproot.open("../data/TTbar_pu200_D49_prompt.root")["L1TrackNtuple/eventTree"].arrays("*", namedecode="utf-8"))

# fake rate plots
bins, bins_err, faker_disp, faker_err_disp = faker_plot(arrays_disp, "pt")
bins, bins_err, faker_pro, faker_err_pro = faker_plot(arrays_pro, "pt")
plt.style.use(hep.style.CMS)
plt.errorbar(bins,faker_disp,xerr=bins_err,yerr=faker_err_disp,fmt='o',label="ttbar displaced")
plt.errorbar(bins,faker_pro,xerr=bins_err,yerr=faker_err_pro,fmt='o',label="ttbar prompt")
plt.legend(loc="best")
plt.xlabel(r"L1 track $p_T$ (GeV/c)")
plt.ylabel("fake rate")
plt.ylim(0,1.1)
plt.show()

bins, bins_err, faker_disp, faker_err_disp = faker_plot(arrays_disp, "eta")
bins, bins_err, faker_pro, faker_err_pro = faker_plot(arrays_pro, "eta")
plt.style.use(hep.style.CMS)
plt.errorbar(bins,faker_disp,xerr=bins_err,yerr=faker_err_disp,fmt='o',label="ttbar displaced")
plt.errorbar(bins,faker_pro,xerr=bins_err,yerr=faker_err_pro,fmt='o',label="ttbar prompt")
plt.legend(loc="best")
plt.xlabel(r"L1 track $\eta$")
plt.ylabel("fake rate")
plt.ylim(0,1.1)
plt.show()

bins, bins_err, faker_disp, faker_err_disp = faker_plot(arrays_disp, "phi")
bins, bins_err, faker_pro, faker_err_pro = faker_plot(arrays_pro, "phi")
plt.style.use(hep.style.CMS)
plt.errorbar(bins,faker_disp,xerr=bins_err,yerr=faker_err_disp,fmt='o',label="ttbar displaced")
plt.errorbar(bins,faker_pro,xerr=bins_err,yerr=faker_err_pro,fmt='o',label="ttbar prompt")
plt.legend(loc="best")
plt.xlabel(r"L1 track $\phi$")
plt.ylabel("fake rate")
plt.ylim(0,1.1)
plt.show()

bins, bins_err, faker_disp, faker_err_disp = faker_plot(arrays_disp, "z0")
bins, bins_err, faker_pro, faker_err_pro = faker_plot(arrays_pro, "z0")
plt.style.use(hep.style.CMS)
plt.errorbar(bins,faker_disp,xerr=bins_err,yerr=faker_err_disp,fmt='o',label="ttbar displaced")
plt.errorbar(bins,faker_pro,xerr=bins_err,yerr=faker_err_pro,fmt='o',label="ttbar prompt")
plt.legend(loc="best")
plt.xlabel(r"L1 track $z_0$ (cm)")
plt.ylabel("fake rate")
plt.ylim(0,1.1)
plt.show()

bins, bins_err, faker_disp, faker_err_disp = faker_plot(arrays_disp, "d0")
#bins, bins_err, faker_pro, faker_err_pro = faker_plot(arrays_pro, "d0")
plt.style.use(hep.style.CMS)
plt.errorbar(bins,faker_disp,xerr=bins_err,yerr=faker_err_disp,fmt='o',label="ttbar displaced")
#plt.errorbar(bins,faker_pro,xerr=bins_err,yerr=faker_err_pro,fmt='o',label="ttbar prompt")
plt.legend(loc="best")
plt.xlabel(r"L1 track $|d_0|$ (cm)")
plt.ylabel("fake rate")
plt.ylim(0,1.1)
plt.show()