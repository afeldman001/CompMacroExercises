# ================================================== #
# This script provides a replication of the baseline # 
# RBC model with basic government spending, as shown # 
# in Baxter & King (1993, AER), Section III. The     #
# program offers computation of steady state values, #
# simulates the macroeconomic impact of permanent    #
# government purchases (fiscal shock),and reproduces #
# Figure II from the paper.                          #
# ================================================== #

# import libraries
import sympy as sp
import numpy as np
import pandas as pd
from writeOut import writeOut

### Declare symbolic endogenous variables ###

# model name
model_name = "rbc"

# declare endogenous variables
endo_names = [
    "y", "c", "l", "n", "iv", "k", "lam", "tr", "tau", "gb",
    "w", "q", "r", "uc", "ul", "fn", "fk", "check_walras"
]

# declare exogenous variables
exo_names = ["e_gb"]

# declare parameters 
param_names = [
    "A", "GAMMAX", "BETA", "DELTA_K", "THETA_L", 
    "THETA_K", "THETA_N", "GB_BAR", "TAU_BAR"
]

# compute lengths for later use
endo_nbr = len(endo_names)
exo_nbr = len(exo_names)
param_nbr = len(param_names)

# time indexes
time_ind = ["tm1", "t", "tp1", "ss"]

# map time indexes to symbolic endogenous variables
symbols = {}
for var in endo_names:
    for ind in time_ind:
        symbols[f"{var}_{ind}"] = sp.Symbol(f"{var}_{ind}")

# exogenous variables
for eps in exo_names:
    symbols[eps] = sp.Symbol(eps)

# parameters 
for param in param_names:
    symbols[param] = sp.Symbol(param)

### Model equations using symbolic mapping ###

# declare shorthand 
s = symbols 

# dynamic model equations
dynamic_eqs = [
    s["uc_t"] - s["c_t"]**(-1),
    s["ul_t"] - s["THETA_L"] * s["l_t"]**(-1),
    s["y_t"] - s["A"] * s["k_tm1"]**s["THETA_K"] * s["n_t"]**s["THETA_N"],
    s["fk_t"] - s["THETA_K"] * s["A"] * s["k_tm1"]**(s["THETA_K"] - 1) * s["n_t"]**s["THETA_N"],
    s["fn_t"] - s["THETA_N"] * s["A"] * s["k_tm1"]**s["THETA_K"] * s["n_t"]**(s["THETA_N"] - 1),
    s["GAMMAX"] * s["k_t"] - ((1 - s["DELTA_K"]) * s["k_tm1"] + s["iv_t"]),
    s["l_t"] + s["n_t"] - 1,
    s["c_t"] + s["iv_t"] - (1 - s["tau_t"]) * s["y_t"] - s["tr_t"] - s["check_walras_t"],
    s["c_t"] + s["iv_t"] + s["gb_t"] - s["y_t"],
    s["tau_t"] * s["y_t"] - s["gb_t"] - s["tr_t"],
    s["uc_t"] - s["lam_t"],
    s["ul_t"] - s["lam_t"] * (1 - s["tau_t"]) * s["fn_t"],
    s["BETA"] * s["lam_tp1"] * (s["q_tp1"] + 1 - s["DELTA_K"]) - s["GAMMAX"] * s["lam_t"],
    s["q_t"] - (1 - s["tau_t"]) * s["fk_t"],
    s["gb_t"] - (s["GB_BAR"] + s["e_gb"]),
    s["tau_t"] - s["TAU_BAR"],
    1 + s["r_t"] - (s["GAMMAX"] * s["lam_t"] / (s["lam_tp1"] * s["BETA"])),
    s["w_t"] - s["fn_t"]
]

if len(dynamic_eqs) != endo_nbr:
    raise ValueError("You need to have as many endogenous variables as model equations. BAD!")

### Create lead-lag incidence matrix for dynamic variables ###

# declare empty matrix
lead_lag_incidence = np.zeros((3, endo_nbr), dtype=int)
idx = 1

# loop through endogenous variable names and checks if time_indexed variables
# appear in dynamic_eqs 
for j, name in enumerate(endo_names):
    if f"{name}_tm1" in symbols and any(symbols[f"{name}_tm1"] in eq.free_symbols for eq in dynamic_eqs):
        lead_lag_incidence[0, j] = idx
        idx += 1
    if f"{name}_t" in symbols and any(symbols[f"{name}_t"] in eq.free_symbols for eq in dynamic_eqs):
        lead_lag_incidence[1, j] = idx
        idx += 1
    if f"{name}_tp1" in symbols and any(symbols[f"{name}_tp1"] in eq.free_symbols for eq in dynamic_eqs):
        lead_lag_incidence[2, j] = idx
        idx += 1

# distinguish endogenous variables by timing structure #
lead = lead_lag_incidence[2, :] != 0
lag  = lead_lag_incidence[0, :] != 0
curr = lead_lag_incidence[1, :] != 0

endo_static_names = [endo_names[i] for i in range(endo_nbr) if not lag[i] and curr[i] and not lead[i]]
endo_pred_names   = [endo_names[i] for i in range(endo_nbr) if     lag[i] and not lead[i]]
endo_fwrd_names   = [endo_names[i] for i in range(endo_nbr) if not lag[i] and     lead[i]]
endo_mixed_names  = [endo_names[i] for i in range(endo_nbr) if     lag[i] and     lead[i]]

# ordered variable names as in Dynare: static, pred, mixed, fwrd
order_var = endo_static_names + endo_pred_names + endo_mixed_names + endo_fwrd_names

# construct ordered dynamic variable list based on incidence matrix
dynamic_names = (
    [f"{endo_names[i]}_tm1" for i in range(endo_nbr) if lead_lag_incidence[0, i] != 0] +
    [f"{endo_names[i]}_t"   for i in range(endo_nbr) if lead_lag_incidence[1, i] != 0] +
    [f"{endo_names[i]}_tp1" for i in range(endo_nbr) if lead_lag_incidence[2, i] != 0]
)
dynamic_syms = [symbols[name] for name in dynamic_names]

# display lead-lag incidence matrix for inspection 
lead_lag_df = pd.DataFrame(
    lead_lag_incidence,
    index=["t-1", "t", "t+1"],
    columns=endo_names
)
print("\nLead-Lag Incidence Matrix:")
print(lead_lag_df)

### Compute dynamic and static Jacobian matrices ###

# compute dynamic Jacobian matrix
exo_syms = [symbols[name] for name in exo_names]
dynamic_eqs = sp.Matrix(dynamic_eqs) # convert to symbolic matrix
dynamic_g1 = sp.Matrix(dynamic_eqs).jacobian(dynamic_syms + exo_syms)

# generate static model equations (remove time indices)
# replace time indexes and ss sufixed variables with static counterparts
static_eqs = []

for eq in dynamic_eqs:
    e = eq
    for var in endo_names:
        for suffix in ["tm1", "t", "tp1", "ss"]:
            e = e.subs(symbols[f"{var}_{suffix}"], sp.Symbol(var))
    static_eqs.append(e)

static_eqs = sp.Matrix(static_eqs)  # create column vector matrix

# list of static (unsuffixed) endogenous symbols
static_syms = [sp.Symbol(name) for name in endo_names]

# compute static Jacobian
static_g1 = sp.Matrix(static_eqs).jacobian(static_syms)

### Write symbolic expressions and Jacobians to Python script files ###

writeOut(static_eqs, model_name + "_static_resid", 'residual', is_static=True,
         dynamic_names=dynamic_names, endo_names=endo_names,
         exo_names=exo_names, param_names=param_names)

writeOut(static_g1, model_name + "_static_g1", 'g1', is_static=True,
         dynamic_names=dynamic_names, endo_names=endo_names,
         exo_names=exo_names, param_names=param_names)

writeOut(dynamic_eqs, model_name + "_dynamic_resid", 'residual', is_static=False,
         dynamic_names=dynamic_names, endo_names=endo_names,
         exo_names=exo_names, param_names=param_names)

writeOut(dynamic_g1, model_name + "_dynamic_g1", 'g1', is_static=False,
         dynamic_names=dynamic_names, endo_names=endo_names,
         exo_names=exo_names, param_names=param_names)

# store to structure whcih can be passed to solvers,
# or used to auto-generate simulation routines
MODEL = {
    "fname": model_name,
    "endo_names": endo_names,
    "endo_nbr": endo_nbr,
    "nstatic": len(endo_static_names),
    "npred": len(endo_pred_names),
    "nboth": len(endo_mixed_names),
    "nfwrd": len(endo_fwrd_names),
    "nspred": len(endo_pred_names) + len(endo_mixed_names),
    "nsfwrd": len(endo_mixed_names) + len(endo_fwrd_names),
    "exo_names": exo_names,
    "exo_nbr": exo_nbr,
    "lead_lag_incidence": lead_lag_incidence,
    "param_names": param_names,
    "param_nbr": param_nbr,
    "order_var": order_var
}
