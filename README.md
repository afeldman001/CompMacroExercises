# CompMacroExercises

This repo is for R-based replications of macroeconomic exercises from my Computational Macroeconomics course. Most of the original problems are in Dynare or MATLAB, I’m translating them into R to better understand the models and solvers.

The current version of the perfect foresight solver implements a full-system Newton-Raphson method using dense numerical Jacobians and stacked residuals. It reproduces the basic logic of Dynare’s solver but omits symbolic derivation, sparse linear algebra, and homotopy continuation. Future revisions will address these limitations.
