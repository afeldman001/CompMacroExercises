# ---------------------------------------------------------------------------- #

# perfect_foresight_solver.R
# modular perfect foresight path solver for DSGE models in R
# compatible with time-indexed Baxter & King (1993)-style DSGE models

# Simulate transition using perfect foresight: solve nonlinear system over time
# NOTE: This implementation uses a custom Newton-style solver that stacks the full system
# across time periods and solves it as a nonlinear root-finding problem.
# While it captures the general logic of Dynare's `perfect_foresight_solver`, it lacks:
# - Sparse Jacobian support (for efficiency on large models)
# - Homotopy continuation (for improving convergence on hard nonlinearities)
# - Advanced fallback options (e.g., linear or marginal linearization)
# As such, this version is slower and less robust than Dynare’s implementation.
# The purpose here is expository: to gain hands-on practice with numerical methods
# and understand the structure of deterministic simulations in DSGE models.

# ---------------------------------------------------------------------------- #

perfect_foresight_solver <- function(
    model_equations,        # function(vars_t, vars_tm1, vars_tp1, params, exo_t)
    init_path,              # matrix: T x n_vars, initial guess for endogenous vars
    params,                 # list of parameters
    exo_path,               # matrix: T x n_exo, values for exogenous vars
    var_names,              # character vector: names of endogenous vars
    exo_names,              # character vector: names of exogenous vars
    steady_state_init,      # named vector for period 0 (t=0) steady state
    steady_state_terminal,  # named vector for period T+1 (t=T+1) steady state
    method = "Newton",      # solver method
    tol = 1e-8,             # tolerance
    max_iter = 500          # max iterations
) {
  T <- nrow(init_path)
  n_vars <- length(var_names)
  
  # flatten initial guess matrix into vector
  x_vec <- as.vector(t(init_path))
  
  # residual function for entire time path
  residual_fn <- function(x) {
    X <- matrix(x, nrow = T, byrow = TRUE)
    res <- numeric(T * n_vars)
    
    for (t in 1:T) {
      vars_t   <- setNames(X[t, ], var_names)
      
      vars_tm1 <- setNames(
        if (t == 1) steady_state_init[var_names] else X[t - 1, ],
        paste0(var_names, "_tm1")
      )
      
      vars_tp1 <- setNames(
        if (t == T) steady_state_terminal[var_names] else X[t + 1, ],
        paste0(var_names, "_tp1")
      )
      
      exo_t    <- setNames(exo_path[t, ], exo_names)
      
      res[((t - 1) * n_vars + 1):(t * n_vars)] <-
        model_equations(vars_t, vars_tm1, vars_tp1, params, exo_t)
    }
    
    return(res)
  }
  
  start_time <- Sys.time()
  
  # solve the stacked nonlinear system
  result <- nleqslv::nleqslv(
    x = x_vec,
    fn = residual_fn,
    method = method,
    control = list(ftol = tol, maxit = max_iter)
  )
  
  end_time <- Sys.time()
  elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  # return path in matrix form
  X_solution <- matrix(result$x, nrow = T, byrow = TRUE)
  colnames(X_solution) <- var_names
  
  # compute residuals for diagnostics
  final_residuals <- residual_fn(result$x)
  residual_norm <- max(abs(final_residuals))
  
  # print diagnostics (optional)
  cat("\n--- Perfect Foresight Solver Summary ---\n")
  cat("Converged:       ", result$termcd == 1, "\n")
  cat("Iterations:      ", result$iter, "\n")
  cat("Max residual:    ", format(residual_norm, digits = 3), "\n")
  cat("Elapsed time:    ", round(elapsed_time, 3), "seconds\n")
  cat("----------------------------------------\n")
  
  # return full output
  list(
    path = X_solution,
    converged = result$termcd == 1,
    iterations = result$iter,
    residual_norm = residual_norm,
    residuals = final_residuals,
    time_seconds = elapsed_time,
    solver_output = result
  )
}
