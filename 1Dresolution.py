import uproot
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep

# take in .root into arrays
arrays_disp = (uproot.open("../data/TTbar_pu200_D49_extended.root")["L1TrackNtuple/eventTree"].arrays("*", namedecode="utf-8"))
arrays_pro = (uproot.open("../data/TTbar_pu200_D49_prompt.root")["L1TrackNtuple/eventTree"].arrays("*", namedecode="utf-8"))

# 1D resolution plots
plt.style.use(hep.style.CMS)
plt.hist(arrays_disp["tp_pt"].flatten()-arrays_disp["matchtrk_pt"].flatten(),bins=50,range=(-50,50),histtype="step",label="ttbar displaced")
plt.hist(arrays_pro["tp_pt"].flatten()-arrays_pro["matchtrk_pt"].flatten(),bins=50,range=(-50,50),histtype="step",label="ttbar prompt")
plt.legend(loc="best")
plt.xlabel(r"$p_T$ resolution (GeV/c)")
plt.ylabel("counts")
plt.yscale("log")
plt.show()

plt.style.use(hep.style.CMS)
plt.hist(arrays_disp["tp_eta"].flatten()-arrays_disp["matchtrk_eta"].flatten(),bins=50,range=(-.5,.5),histtype="step",label="ttbar displaced")
plt.hist(arrays_pro["tp_eta"].flatten()-arrays_pro["matchtrk_eta"].flatten(),bins=50,range=(-.5,.5),histtype="step",label="ttbar prompt")
plt.legend(loc="best")
plt.xlabel(r"$\eta$ resolution")
plt.ylabel("counts")
plt.yscale("log")
plt.show()

plt.style.use(hep.style.CMS)
plt.hist(arrays_disp["tp_phi"].flatten()-arrays_disp["matchtrk_phi"].flatten(),bins=50,range=(-.5,.5),histtype="step",label="ttbar displaced")
plt.hist(arrays_pro["tp_phi"].flatten()-arrays_pro["matchtrk_phi"].flatten(),bins=50,range=(-.5,.5),histtype="step",label="ttbar prompt")
plt.legend(loc="best")
plt.xlabel(r"$\phi$ resolution")
plt.ylabel("counts")
plt.yscale("log")
plt.show()

plt.style.use(hep.style.CMS)
plt.hist(arrays_disp["tp_z0"].flatten()-arrays_disp["matchtrk_z0"].flatten(),bins=50,range=(-.5,.5),histtype="step",label="ttbar displaced")
plt.hist(arrays_pro["tp_z0"].flatten()-arrays_pro["matchtrk_z0"].flatten(),bins=50,range=(-.5,.5),histtype="step",label="ttbar prompt")
plt.legend(loc="best")
plt.xlabel(r"$z_0$ resolution (cm)")
plt.ylabel("counts")
plt.yscale("log")
plt.show()

plt.style.use(hep.style.CMS)
plt.hist(arrays_disp["tp_d0"].flatten()-arrays_disp["matchtrk_d0"].flatten(),bins=50,range=(-10,10),histtype="step",label="ttbar displaced")
#plt.hist(arrays_pro["tp_d0"].flatten()-arrays_pro["matchtrk_d0"].flatten(),bins=50,range=(-10,10),histtype="step",label="ttbar prompt")
plt.legend(loc="best")
plt.xlabel(r"$d_0$ resolution (cm)")
plt.ylabel("counts")
plt.yscale("log")
plt.show()