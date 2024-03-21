"""
Manually restart the line search
"""
optimize._line_search._restart_line_search()
optimize.initialize_search()
alpha, _ = optimize.calculate_step_length()

alpha = 0.5

m_try = optimize.compute_trial_model(alpha=alpha)

optimize.save_vector(name="m_try", m=m_try)
optimize.save_vector(name="alpha", m=alpha)
# optimize.checkpoint()

# Expose model `m_try` to the solver by placing it in eval_func dir.
_path_m_try = os.path.join(workflow.path.eval_func, "model")
m_try.write(path=_path_m_try)
logger.info(f"`m_try` model parameters for initial line search step")
solver.check_model_values(path=_path_m_try)

system.array = 0
workflow.evaluate_line_search_misfit()
