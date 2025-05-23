# ============================================================================
# compute_steady_state.R
# Dynare-style steady state solver for DSGE models in R
# ============================================================================
# This function takes:
# - A user-defined steady-state residual function (model_equations_ss),
# - An initial guess for the endogenous variables (initval),
# - Parameter list and exogenous values,
# And returns:
# - A named vector of steady-state endogenous values
# - A convergence flag and full solver diagnostics
# ============================================================================

compute_steady_state <- function(
    model_equations_ss,  # function(vars, params, exo): returns residuals
    init_guess,          # named numeric vector: initial values of endogenous vars
    params,              # named list of parameters
    exo_vals,            # named list or vector of exogenous values (e.g., e_gb = 0)
    method = "Broyden",  # nleqslv method (e.g., Broyden, Newton)
    tol = 1e-10,         # solver tolerance
    max_iter = 500       # max solver iterations
) {
  # residual wrapper for nleqslv
  residual_fn <- function(x) {
    names(x) <- names(init_guess)  # preserve names
    model_equations_ss(x, params, exo_vals)
  }
  
  # solve the steady state system
  result <- nleqslv::nleqslv(
    x = init_guess,
    fn = residual_fn,
    method = method,
    control = list(ftol = tol, maxit = max_iter)
  )
  
  # return results
  list(
    steady_state = setNames(result$x, names(init_guess)),
    converged = result$termcd == 1,
    info = result
  )
}
