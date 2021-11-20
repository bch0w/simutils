"""
Nov. 19 2021 

Pyaflowa multiprocessing implementation for Emanuele Casarottis data,
testing out the parallel processing functionalities of Pyaflowa.

Problem at hand: 
    - Synthetics have already been generated and are sitting in a SPECFEM
    OUTPUT_FILES directory
    - Observed data, event metdata, station metadata are already stored in
    PyASDF ASDFDataSets, data do NOT need to be gathered via FDSN
    - Large event number (N~=200), processing should be done in parallel,
    using Pyaflowas internal calls to concurrent.futures
    - Desired outputs are adjoint sources and 
"""
import os
import sys
import shutil
from pyatoa import Pyaflowa, Config

# SETUP EXAMPLE PROBLEM
if not os.path.exists("./ITALY_TRAINING/772141_1.h5"):
    print("Setting up directory, copying ASDFDataSet to get 'multiple events'")
    shutil.copyfile(src="./ITALY_TRAINING/772141_2.h5",
                    dst="./ITALY_TRAINING/772141_1.h5"
                    )
    shutil.copyfile(src="./ITALY_TRAINING/772141_2.h5",
                    dst="./ITALY_TRAINING/772141_3.h5"
                    )
    print("Symlinking synthetic data to get 'multiple events'")
    for src in ["772141_1", "772141_2", "772141_3"]:
        os.symlink(src="./ascii_seismograms",
                   dst=f"./OUTPUT_FILES/{src}")

# !!! PARAMETER SET
mode = "parallel"
max_workers = 4  # Maximum number of parallel processes to run
cfg = Config(iteration=1, step_count=0, min_period=10., max_period=50.,
             filter_corners=2, unit_ouput="VEL", pyflex_preset="nzni1D_10-30s",
             start_pad=100, client=None)

# Tell Pyaflowa where your data is and where outputs should go
paths = {"workdir": "./",   # CWD
         "synthetics": "./OUTPUT_FILES/{source_name}",  # SPECFEM3D format
         "stations_file": "./STATIONS_TRAINING/STATIONS_772141_2",  # STATIONS
         "datasets": "./ITALY_TRAINING"  # ASDFDataSets
         }

# Casarotti data has just been copied a few times to get 'multiple events'
source_names = ["772141_1", "772141_2", "772141_3"]

# !!! PARAMETER SET

# Run in standalone mode to get create a local directory structure
pf = Pyaflowa(structure="standalone", config=cfg, log_level="DEBUG", **paths)

if mode == "serial":
    for source_name in source_names:
        misfit = pf.process_event(source_name)
elif mode == "parallel":
    misfits = pf.multi_event_process(source_names=source_names,
                                     max_workers=max_workers
                                     )


