import matplotlib.pyplot as plt
from pyasdf import ASDFDataSet
from pyatoa import Manager
from pyflex import logger
from IPython import embed

logger.setLevel("DEBUG")


iter_ = "i07/s00"
min_periods = [20]
max_period = 30
sta_check = "AK.COLD"
tshift = 0
sta_check = None

ds = ASDFDataSet("AK_A21K.h5")
mgmt = Manager(ds=ds)
nwin = 0

for min_period in min_periods:
    for sta in ds.waveforms.list():
        if sta_check is not None:
            if sta != sta_check:
                continue

        try:
            print(sta, end="... ")

            mgmt.load(sta, iter_)

            mgmt.config.min_period = min_period
            mgmt.config.max_period = max_period
            mgmt.config.component_list = ["Z"]

            mgmt.config.pyflex_config.cc_acceptance_level = 0.6
            mgmt.config.pyflex_config.tshift_acceptance_level = 10
            mgmt.config.pyflex_config.tshift_reference = 0

            mgmt.st_obs = mgmt.st_obs.select(component="Z")
            mgmt.st_syn = mgmt.st_syn.select(component="Z")

            mgmt.standardize()
            mgmt.preprocess(normalize_to="syn")

            if False:
                plt.plot(mgmt.st_obs[0].times(), mgmt.st_obs[0].data, "k")
                plt.plot(mgmt.st_syn[0].times() + tshift, mgmt.st_syn[0].data, "r")
                plt.xlim([75, 300])
                plt.savefig(f"figures/custom_{sta}_{min_period}_{max_period}.png")
                plt.close("all")

            mgmt.window()
            mgmt.measure()

            mgmt.plot(save=f"figures/{sta}_{min_period}_{max_period}.png", show=False)
            plt.close("all")

            nwin += mgmt.stats.nwin

            if sta_check is not None:
                break
            print("")
        except Exception as e:
            print(e)
            continue

    print(f"\nWINDOWS = {nwin}")

