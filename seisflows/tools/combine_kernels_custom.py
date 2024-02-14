"""
Custom SeisFlows Workflow functions to run tasks on system.
Copy paste into your Workflow source code and run function from debugger.
"""
    def custom_combine_kernels(self):
        """
        Generate event RR and TT kernels by summing individual E? and N?
        kernels for each source. 
        """
        def combine_kernels(**kwargs):
            parameters = [f"{par}_kernel" for par in self.solver._parameters]
            for kernel in ["R", "T"]:
                for srcname in self.solver.source_names:
                    base_path = os.path.join(self.path.eval_grad, 
                                             "kernels", srcname)

                    input_paths = [os.path.join(base_path, f"E{kernel}"),
                                   os.path.join(base_path, f"N{kernel}")]
                    output_path = os.path.join(base_path, f"{kernel}{kernel}")

                    self.solver.combine(input_paths, output_path, parameters)
        self.system.run([combine_kernels], single=True)

    def custom_function_sum_kernels(self):
        """
        Generate ZZ, RR or TT misfit kernels through summation of misfit kernels
        """
        def sum_kernels(**kwargs):
            kernels = ["ZZ"]  # For after 'evaluate_zz_adjoint_simulation'
            # kernels = ["ZZ", "RR", "TT"]
            parameters = [f"{par}_kernel" for par in self.solver._parameters]
            base_path = os.path.join(self.path.eval_grad, "kernels")
            for kernel in kernels:
                input_paths = []
                output_path = os.path.join(base_path, f"{kernel}")

                for srcname in self.solver.source_names:
                    input_paths.append(os.path.join(base_path, srcname, kernel))

                self.solver.combine(input_paths, output_path, parameters)
        self.system.run([sum_kernels], single=True)

