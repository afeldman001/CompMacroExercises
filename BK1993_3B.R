# ============================================================================ #
#                                                                              #
# This script provides a replication of the Baxter & King (1993, AER)          #
# Benchmark Model with Basic Government Purchases Figure 2 from section        #
# IIIB. This script ports Willi. Mutschler's (2023, University of              #  
# Tuebingen) version in R (I add the term structure) as a practice             #
# exercise.                                                                    #
# ---------------------------------------------------------------------------- #                                                                                                                                     #
# Aaron Feldman                                                                #
# April 30, 2025                                                               #
# University of Tübingen                                                       #
# Computational Macroeconomics SS25                                            #
#                                                                              #
# ============================================================================ #

##### Housekeeping #####

# clear environment
rm(list = ls())

# clear console
cat("\014")  # sends form feed to console to simulate clearing

##### Install Required Packages and Load Libraries #####

# install required packages if not already installed
packages <- c("tidyverse", "nleqslv", "Matrix", "numDeriv")
installed <- packages %in% installed.packages()
if (any(!installed)) install.packages(packages[!installed])

# load libraries
library(tidyverse)  # data manipulation and plotting (includes dplyr, ggplot2, etc.)
library(nleqslv)    # solving systems of nonlinear equations (for steady-state and transition paths) 
library(Matrix)     # matrix algebra and numerical routines
library(numDeriv)   # for Jacobian computation  

source("perfect_foresight_solver.R")
source("compute_steady_state.R")

##### 1. Declare Variables and Parameters #####

# endogenous variables 

variables <- c(
  "y",       # output
  "c",       # consumption
  "l",       # leisure
  "n",       # labor
  "iv",      # private investment
  "k",       # private capital stock
  "lam",     # Lagrangian multiplier on budget constraint (equation 7 in BK (1993))
  "tr",      # fiscal transfers
  "tau",     # net tax rate
  "gb",      # basic government spending
  "w",       # real wage
  "q",       # after tax rental rate of capital
  "r",       # real interest rate
  "uc",      # marginal utility wrt consumption
  "ul",      # marginal utility wrt leisure
  "fn",      # marginal product of labor
  "fk",      # marginal product of capital
  "check_walras" # redundant resource constraint checker
)

# exogenous variable (shock term)

exogenous <- c(
  "e_gb"  # basic government spending shock
)

# parameters

parameters <- list(
  A = NA,        # common technology level
  GAMMAX = NA,   # growth rate of labor augmenting technology (footnote 3 (BK, 1993))
  BETA = NA,     # discount factor
  DELTA_K = NA,  # capital depreciation rate
  THETA_L = NA,  # weight of leisure in utility function
  THETA_K = NA,  # private capital productivity in production function
  THETA_N = NA,  # labor productivity in production function
  GB_BAR = NA,   # target value of basic government spending
  TAU_BAR = NA  # target tax rate
)

##### 2. Set Model Equations (residuals at time t) #####

model_equations <- function(vars_t, vars_tm1, vars_tp1, params, exo_t) {
  with(c(vars_t, vars_tm1, vars_tp1, params, exo_t), {
    
    eqs <- numeric(18)   # 18 model equations, since 18 endogenous variables
    
    # 1. marginal utility wrt consumption given momentary utility in eq . 2 
    eqs[1] <- uc - c^(-1) # uc = c^(-1)
    
    # 2. marginal utility wrt to leisure given momentary utility in eq. 2
    eqs[2] <- ul - THETA_L * l^(-1) # ul = THETA_L * l^(-1)
    
    # 3. production function (eq. 3)
    eqs[3] <- y - A * k_tm1^THETA_K * n^THETA_N # y = A * k(-1)^THETA_K * n^THETA_N
    
    # 4. marginal product of capital given Cobb-Douglas production function in eq. 3
    eqs[4] <- fk - THETA_K * A * k_tm1^(THETA_K - 1) * n^THETA_N # fk = THETA_K * A * k(-1)^(THETA_K - 1) * n^THETA_N
    
    # 5. marginal product of labor given Cobb-Douglas production function in eq. 3
    eqs[5] <- fn - THETA_N * A * k_tm1^THETA_K * n^(THETA_N - 1) # fn = THETA_N * A * k(-1)^THETA_K * n^(THETA_N - 1)
    
    # 6. law of motion of private capital stock, eq. 4 combined with footnote 3
    eqs[6] <- GAMMAX * k - ((1 - DELTA_K) * k_tm1 + iv) # GAMMAX * k = (1 - DELTA_K) * k(-1) + iv
    
    # 7. time endowment, eq. 6
    eqs[7] <- l + n - 1 # l + n = 1
    
    # 8. resource constraint, eq. 7 (redundant due to Walra's Law)
    eqs[8] <- c + iv - ((1 - tau) * y + tr + check_walras) # c + iv = (1 - tau) * y + tr + check_walras
    
    # 9. market clearing, eq. 8
    eqs[9] <- c + iv + gb - y # c + iv + gb = y
    
    # 10. government budget constraint, eq. 9
    eqs[10] <- tau * y - gb + tr # tau * y = gb + tr
    
    # 11. shadow price of consumption utility, eq. 10
    eqs[11] <- uc - lam # uc = lam
    
    # 12. labor-leisure decision combined with labor demand, eq. 11
    eqs[12] <- ul - lam * (1 - tau) * fn # ul = lam * (1 - tau) * fn
    
    # 13. savings decision, eq. 12 combined with footnote 3
    eqs[13] <- BETA * lam_tp1 * (q_tp1 + 1 - DELTA_K) - GAMMAX * lam 
    
    # 14. after tax rental rate of capital, below eq. 13
    eqs[14] <- q - (1 - tau) * fk 
    
    # 15. fiscal rule: government spending
    eqs[15] <- gb - (GB_BAR + e_gb)
    
    # 16. fiscal rule: tax rate
    eqs[16] <- tau - TAU_BAR
    
    # 17. one period real interest rate (see page 322)
    eqs[17] <- 1 + r - (GAMMAX * lam / (lam_tp1 * BETA)) # 1 + r = GAMMAX * lam / (lam(+1) * BETA)
    
    # 18. real wage = marginal product of labor ("King, Plosser, and Robelo (1988) page 204 above section 3")
    eqs[18] <- w - fn 
    
    return(eqs)
    
  })
}

##### 3. Parameter Calibration (annual) based on benchmark calibration from BK, 1993 (Table 1) #####

# structural parameters 
params <- list()

# normalized technology level
params$A <- 1  

# labor-augmenting productivity growth rate (footnote 3), follows King, Plosser, and Rebello (1988) or King and Rebello (1999)
params$GAMMAX <- 1 + 0.016 

# production function parameters
params$THETA_K <- 0.42 # capital share
params$THETA_N <- 1 - params$THETA_K # labor share (constant returns to scale)
params$DELTA_K <- 0.10 # depreciation rate

# preferences
N <- 0.20   # steady-state labor supply 
L <- 1 - N  # steady-state leisure
R <- 0.065  # steady-state real interest rate 
params$BETA <- params$GAMMAX / (1 + R)  # implied steady-state discount factor

# government
sG <- 0.20              # G/Y ration (government spending share of output)
params$TAU_BAR <- 0.20  # same tax rate as government spending implies zero transfers
TR <- 0                 # zero transfers is the benchmark
TAU <- params$TAU_BAR

# calibrate GB_BAR and THETA_L implicitly via steady-state relationships
Q  <- params$GAMMAX / params$BETA - 1 + params$DELTA_K      # implied after-tax return on capital
FK <- Q / (1 - TAU)                                         # marginal product of capital

K <- (FK / (params$THETA_K * params$A * N^params$THETA_N))^(1 / (params$THETA_K - 1))  # capital demand

FN <- params$THETA_N * params$A * K^params$THETA_K * N^(params$THETA_N - 1)  # marginal product of labor
W  <- FN                                                                     # real wage
IV <- (params$GAMMAX - 1 + params$DELTA_K) * K                               # investment
Y  <- params$A * N^(1 - params$THETA_K) * K^params$THETA_K                   # output
params$GB_BAR <- sG * Y                                                      # target gov spending
GB <- params$GB_BAR                                                          # actual gov spending
C  <- Y - IV - GB                                                            # consumption
UC <- C^(-1)                                                                 # marginal utility of consumption
LAM <- UC                                                                    # shadow price of utility
UL <- LAM * (1 - TAU) * FN                                                   # marginal utility of leisure
params$THETA_L <- UL * L                                                     # implied leisure weight in utility


##### 4. Simulation: Transition to New Steady-State #####

# set initial conditions 
initval <- list(
  y = Y,
  c = C,
  l = L,
  n = N,
  iv = IV,
  k = K,
  lam = LAM,
  tr = TR,
  tau = TAU,
  gb = GB,
  w = W,
  q = Q,
  r = R,
  uc = UC,
  ul = UL,
  fn = FN,
  fk = FK,
  check_walras = 0
)

# wrap existing model equation for steady state (leads/lags set to current)
model_equations_ss <- function(vars, params, exo) {
  vars_tm1 <- setNames(vars, paste0(names(vars), "_tm1")) # use same values for lags
  vars_tp1 <- setNames(vars, paste0(names(vars), "_tp1")) # use same values for leads
  model_equations(
    vars_t = vars,
    vars_tm1 = vars_tm1,
    vars_tp1 = vars_tp1,
    params = params,
    exo_t = exo
  )
}


# compute initial steady state (pre-shock)
init_guess <- unlist(initval) # convert named list to flat vector for solver
exo_vals <- list(e_gb = 0) # no shock in baseline
ss_init <- compute_steady_state(
  model_equations_ss = model_equations_ss,
  init_guess = init_guess,
  params = params,
  exo_vals = exo_vals
)

# set permanent shock: government spending rises by 1% of steady-state output
steady0 <- ss_init$steady_state
commodity_unit <- 0.01 * steady0["y"]     # "commodity unit" (BK, 1993, p. 321)
exo_shock <- list(e_gb = commodity_unit)  # define permanent shock

# recompute new steady state after permanent shock
ss_terminal <- compute_steady_state(
  model_equations_ss = model_equations_ss,
  init_guess = steady0, # use previous SS as guess
  params = params,
  exo_vals = exo_shock
) 

ss_terminal <- setNames(ss_terminal$steady_state, names(steady0))

# set simulation horizon and prepare transition path
T <- 200

# initial guess for the transition path: assume system stays at initial steady state
init_path <- matrix(rep(unlist(steady0), times = T), nrow = T, byrow = TRUE)
colnames(init_path) <- names(steady0)

# set exogenous shock path: permanent increase in government spending from t = 2 onward
exo_path <- matrix(0, nrow = T, ncol = length(exogenous))
colnames(exo_path) <- exogenous
exo_path[1:T, "e_gb"] <- commodity_unit

# simulate transition using perfect foresight (Newton method).
# note: this custom solver lacks Dynare features like sparse Jacobians and homotopy,
# so it's slower and less robust—used here purely for expository purposes.
# warning: convergence takes approximately 4 minutes
solution <- perfect_foresight_solver(
  model_equations = model_equations,
  init_path = init_path,
  params = params,
  exo_path = exo_path,
  var_names = variables,
  exo_names = exogenous,
  steady_state_init = steady0,
  steady_state_terminal = ss_terminal
)

##### 5. Figure 2: Macroeconomic Effects of a Permanent Increase in Basic Government Purchases #####

# extract paths
y_path <- solution$path[, "y"]
c_path <- solution$path[, "c"]
iv_path <- solution$path[, "iv"]
n_path <- solution$path[, "n"]
w_path <- solution$path[, "w"]
r_path <- solution$path[, "r"]
lam_path <- solution$path[, "lam"]

# use steady state from before the shock as baseline
y_ss  <- steady0["y"]
c_ss  <- steady0["c"]
iv_ss <- steady0["iv"]
n_ss  <- steady0["n"]
w_ss  <- steady0["w"]
r_ss  <- steady0["r"]

# compute deviations from steady state
horizon <- 0:20
y_plot  <- (y_path[horizon + 1] - y_ss) / commodity_unit
c_plot  <- (c_path[horizon + 1] - c_ss) / commodity_unit
iv_plot <- (iv_path[horizon + 1] - iv_ss) / commodity_unit
n_plot  <- 100 * (n_path[horizon + 1] / n_ss - 1)
w_plot  <- 100 * (w_path[horizon + 1] / w_ss - 1)
r_plot  <- 10000 * (r_path[horizon + 1] - r_ss)  # in basis points

# compute effect on the term structure (1-20 years)
term_structure <- sapply(1:20, function(n) {
  beta <- params$BETA
  lam_t <- lam_path[2]
  lam_tn <- lam_path[2 + n]
  bond_price <- beta^n * (lam_tn / lam_t)
  yield_post <- -log(bond_price) / n
  yield_ss <- -log(beta)
  10000 * (yield_post - yield_ss)
})

# Figure 2A: Commodity Market
plot(horizon, y_plot, type = "n", ylim = c(-1.2, 1.2), xlim = c(0, 21),
     xlab = "Years", ylab = "Commodity Units", axes = FALSE,
     main = "Figure 2: Commodity Market")
axis(1, at = seq(0, 20, 2))
axis(2, at = seq(-1.2, 1.2, 0.3))
box()
abline(h = 0, lwd = 2)
points(horizon, y_plot, pch = 19)   # output
points(horizon, c_plot, pch = 17)   # consumption
points(horizon, iv_plot, pch = 4)   # investment
legend("bottomright", legend = c("Output", "Consumption", "Investment"),
       pch = c(19, 17, 4), bty = "n")

# Figure 2B: Labor Market
plot(horizon, n_plot, type = "n", ylim = c(-1.6, 1.6), xlim = c(0, 21),
     xlab = "Years", ylab = "Percent", axes = FALSE,
     main = "Figure 2: Labor Market")
axis(1, at = seq(0, 20, 2))
axis(2, at = seq(-1.6, 1.6, 0.4))
box()
abline(h = 0, lwd = 2)
points(horizon, n_plot, pch = 19)  # labor input
points(horizon, w_plot, pch = 17)  # real wage
legend("bottomright", legend = c("Labor Input", "Real Wage"),
       pch = c(19, 17), bty = "n")

# Figure 2C: Financial Market
plot(horizon, r_plot, type = "n", ylim = c(0, 12), xlim = c(0, 21),
     xlab = "Years / Term in Years", ylab = "Basis Points", axes = FALSE,
     main = "Figure 2: Financial Market")
axis(1, at = seq(0, 20, 2))
axis(2, at = seq(0, 12, 1))
box()
abline(h = 0, lwd = 2)
points(horizon, r_plot, pch = 19)             # real interest rate
points(1:20, term_structure, pch = 17)        # term structure at t=1
legend("topright", legend = c("Real Interest Rate", "Term Structure"),
       pch = c(19, 17), bty = "n")


