"""
For debugging and testing updates to preprocessing in an existing SeisFlows
run directory. Copy-paste this into debug mode
"""
# Set parameters here
iteration = 2
step_count = 0
source_name = solver.source_names[0]
component_list = ["R", "T"]
waveform_idx = 3

# Run through setup
preprocess.setup()
obs, syn = preprocess._setup_quantify_misfit(solver.source_name, 
                                             save_adjsrcs=None, 
                                             components="".join(component_list))
# Selects single waveforms to process
o = obs[waveform_idx]
s = syn[waveform_idx]

# Manually run through some setup for the Confib oject
config = preprocess._config.copy()
config.iteration = iteration
config.step_count = step_count
config.event_id = source_name
config.component_list = component_list

# Run preprocessing on login node. Note log will be output in 
# `preprocess.path._tmplogs`
preprocess._quantify_misfit_single(obs_fid=o, syn_fid=s, config=config)

