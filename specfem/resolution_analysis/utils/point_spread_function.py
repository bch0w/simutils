#!/usr/bin/env python
"""
This is the main class seisflows.workflow.sensitivity_analysis
"""
import os
import sys
from glob import glob
import numpy as np
from seisflows.config import custom_import
from seisflows.tools import unix
from seisflows.tools.tools import exists
from seisflows.tools.err import ParameterError
from seisflows.plugins import adjoint

PAR = sys.modules["seisflows_parameters"]
PATH = sys.modules["seisflows_paths"]

system = sys.modules["seisflows_system"]
solver = sys.modules["seisflows_solver"]
preprocess = sys.modules["seisflows_preprocess"]
postprocess = sys.modules["seisflows_postprocess"]


class PointSpreadFunction(custom_import("workflow", "base")):
    """
    Following the methodologies of Fichtner & Trampert 2011b, this class 
    calculates Point Spread Functions using approximate Hessians derived
    by a finite difference calculation of gradients. These kernels can be used
    as a form of resolution testing during adjoint tomography inversions.
    """
    def check(self):
        """ 
        Checks parameters and paths
        """
        # Global defines the location to search for and save scratch files
        if "GLOBAL" not in PATH:
            raise ParameterError(PATH, "GLOBAL")

        if "OUTPUT" not in PATH:
            raise ParameterError(PATH, "OUTPUT")

        # check input
        if "DATA" not in PATH:
            setattr(PATH, 'DATA', None)

        # Will look here for the velocity model
        if "MODEL" not in PATH:
            raise ParameterError(PATH, "MODEL")

        # Here lies the perturbation vector
        if "PERTURB" not in PATH:
            raise ParameterError(PATH, "PERTURB")

        # Initial and True model are not needed, and should point the Model
        if "MODEL_INIT" not in PATH:
            setattr(PATH, 'MODEL_INIT', PATH.MODEL)

        if PATH.DATA is not None and not exists(PATH.DATA):
            assert "MODEL_TRUE" in PATH, "Either 'DATA' or 'MODEL_TRUE'"
            assert exists(PATH.MODEL_TRUE)


    def main(self):
        """ 
        Compute a point spread from a given perturbation 'dm' with respect to a
        reference model 'm'
        """
        # Create the relevant directory structure
        print("SeisFlows Sensitivity Analysis")
        print("\tRunning setup...")
        # self.setup()

        # Load the model and user-provided perturbation
        m = solver.merge(solver.load(PATH.MODEL))
        dm = solver.merge(solver.load(PATH.PERTURB))

        print(f"\tMODEL: {PATH.MODEL}")
        print(f"\tPERTURB: {PATH.PERTURB}")
        print(f"\t\tdm_min={dm.min()}, dm_max={dm.max()}")

        # Set paths to save the gradient caluclations
        path1 = os.path.join(PATH.GLOBAL, "eval1")
        path2 = os.path.join(PATH.GLOBAL, "eval2")

        if not os.path.exists(path1):
            print("\tEvaluating gradient m")
            self.evaluate_gradient(model=m, path=path1)
            unix.cp(PATH.PREPROCESS, os.path.join(path1, "preprocess"))
        else:
            print("\tGradient m already evaluated")

        print("\tEvaluating gradient m + dm")
        # A little hacky but 
        # 1) increment the iteration to 2
        # 2) ensure PAR.FIX_WINDOWS is True 
        # 3) copy-paste datasets from eval1 
        # This will re-evaluate windows from eval1 rather than choosing new ones
        optimize = sys.modules["seisflows_optimize"]
        optimize.iter = 2   
        self.evaluate_gradient(model=(m + dm), path=path2)
        unix.cp(PATH.PREPROCESS, os.path.join(path2, "preprocess"))

        # Retrieve gradients
        path1 = os.path.join(PATH.GLOBAL, "eval1", "gradient")
        path2 = os.path.join(PATH.GLOBAL, "eval2", "gradient")

        grad1 = solver.merge(solver.load(path1, suffix="_kernel"))
        grad2 = solver.merge(solver.load(path2, suffix="_kernel"))

        # Caclulate the action Hessian as the difference of gradients 
        print("\tExporting action Hessian as g(m+dm) - g(m)")
        filename = os.path.join(PATH.OUTPUT, "action_hessian")
        solver.save(save_dict=solver.split(grad2 - grad1), path=filename)


    def setup(self, path=None):
        """
        One-time setup to initiate the SeisFlows working directory
        """
        path = path or PATH.GLOBAL

        # Prepare directory structure
        unix.rm(path)
        unix.mkdir(path)

        # Set up workflow machinery, only need pre and post processing
        preprocess.setup()
        postprocess.setup()

        # Initialize solver directories, generate synthetics etc.
        solver.initialize_solver_directories()

        # Generate 'data' by creating synthetics 
        system.run_single("solver", "setup", model="init")


    def evaluate_gradient(self, model, path):
        """ 
        Performs forward simulation to evaluate objective function and adjoint
        simulation to generate misfit kernels
        """
        unix.mkdir(path)
        solver.save(save_dict=solver.split(model), 
                    path=os.path.join(path, "model"))

        # Evaluate objective function, calculate adjoint sources
        system.run("solver", "eval_fwd", path=path)
        system.run_ancil("solver", "eval_misfit", path=path)

        # Calculate the misfit kernels
        system.run("solver", "eval_grad", path=path)

        # Merge, mask and smooth the gradient
        postprocess.write_gradient(path=path)


