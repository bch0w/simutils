"""
Very simple script to generate a valid MODEL directory from an output SeisFlows
model directory by symlinking the required auxiliary files from MODEL_INIT 
This new output model can then be used as a MODEL_INIT for a new inversion, 
or with the model slicer to create a regular xyz file
"""
import os
import json
from glob import glob


def output_model(output_dir="./output", model_out=None):
    """
    Find MODEL_INIT, symlink all aux files into `model_out`. 
    If model out not defined, will choose the largest model number.

    :type output_dir: str
    :param output_dir: path to the main SeisFlows directory
    :type model_out: int
    :param model_out: model number to be used for output. Defaults to the 
        largest value. Ignores model_true and model_init
    """ 
    # Find the path of MODEL_INIT via the parameter file
    par_file = os.path.join(output_dir, "seisflows_paths.json")
    with open(par_file) as f:
        model_init = json.load(f)["MODEL_INIT"]

    assert(os.path.exists(model_init)), \
                    f"MODEL_INIT does not exist\n{model_init}"
    print(f"MODEL INIT: {model_init}")

    # Determine the model number, only choose numbers, no 'init' or 'true'
    if model_out is None:
        available_models = glob(os.path.join(output_dir, "model_[0-9]???"))
        model_out = sorted(available_models)[-1]
    else:
        model_out = os.path.join(output_dir, model_out)

    assert(os.path.exists(model_out)), f"MODEL_OUT does not exist\n{model_out}"
    print(f"MODEL OUT: {model_out}")

    # Quick check to make sure NPROC is the same for each directory
    nproc_check = [0, 0]
    for i, m in enumerate([model_init, model_out]):
        nprocs = [os.path.basename(_) for _ in glob(os.path.join(m, "*"))]
        # list comprehension strips string parts, e.g. 'proc000001_vp.bin' -> 1
        nproc_check[i] = max([int(_.split('_')[0][4:]) for _ in nprocs])
    assert(nproc_check[0] == nproc_check[1]), f"NPROCS differ {nproc_check}"
    print(f"NPROC: {nproc_check[0]}")
        
    # Symlink all available files that don't already exist in model_out
    model_init_files = glob(os.path.join(model_init, "*"))
    for src in model_init_files:
        dst = os.path.join(model_out, os.path.basename(src))
        if os.path.exists(dst):
            continue
        else:
            os.symlink(src, dst)
        

if __name__ == "__main__":
    output_model()
    
