"""
Utility function for SeisFlows to export trace/waveform data from the SeisFlows
working directory scratch/ subdirectory into a more permanent position on disk.
This is a stop-gap solution for something that should already be in place and
was written to move synthetic Greens function data from an ambient noise adjoint
tomography inversion so that it could be used as target data in a synthetic 
inversion
"""
import os
from glob import glob

input_dir = "/import/c1/ERTHQUAK/bhchow/work/seisflows/GEN_COARSE_DATA/"
output_dir = "/import/c1/ERTHQUAK/bhchow/data/egfs/S40RTS_CRUST1_SGF"
comps = ["R", "T"]
kernels = ["RR", "TT"]

# Start processing here
solver_dir = os.path.join(input_dir, "scratch", "solver")
sources = glob(os.path.join(solver_dir, "*_*"))

for source in sources:
    source_name = os.path.basename(source)
    for comp, kernel in zip(comps, kernels):
        tracedir = os.path.join(source, "traces", "syn", 
                                f"*.MX[{comp}].sem.ascii")
        traces = glob(tracedir)
        for trace in traces:
            trace_name = os.path.basename(trace)
            dst_dir = os.path.join(output_dir, source_name, kernel)
            if not os.path.exists(dst_dir):
                os.makedirs(dst_dir)
            dst = os.path.join(dst_dir, trace_name)
            os.rename(trace, dst)

# input_dir
