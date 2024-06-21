system.partition = "debug"
system.tasktime = 10
system.ntask_max = 1
system.array = 0
preprocess.fix_windows = "ITER"
preprocess.min_period = 25
workflow.setup()
workflow.evaluate_zz_misfit(_preproc_only=True)
