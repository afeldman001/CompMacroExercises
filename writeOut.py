import os
import sympy as sp

def writeOut(Output, name_of_function, name_of_output, is_static,
             dynamic_names, endo_names, exo_names, param_names):
    
    """
    Helper function to write symbolic residuals or Jacobian expressions 
    to a Python function file.
    
    Parameters:
    - Output: SymPy Matrix of symbolic expressions
    - name_of_function: str, base name of the output .py file
    - name_of_output: str, name of the output variable (e.g. 'residual' or 'g1')
    - is_static: bool, whether the model is static (True) or dynamic (False)
    - dynamic_names, endo_names, exo_names, param_names: lists of variable names (str)
    """
    filename = f"{name_of_function}.py"
    
    # delete old version of file (if it exists)
    if os.path.exists(filename):
        os.remove(filename)

    with open(filename, "w") as f:
        # function signature
        if is_static:
            f.write(f"def {name_of_function}(endo_vars, exo_vars, params):\n")
        else:
            f.write(f"def {name_of_function}(dynamic_vars, exo_vars, params, steady_state):\n")

        # indentation helper
        indent = "    "

        # endogenous variable declarations
        if is_static:
            f.write(f"{indent}# Evaluate numerical values for endogenous variables\n")
            for j, name in enumerate(endo_names):
                f.write(f"{indent}{name} = endo_vars[{j}]\n")
        else:
            f.write(f"{indent}# Evaluate numerical values for dynamic variables\n")
            for j, name in enumerate(dynamic_names):
                f.write(f"{indent}{name} = dynamic_vars[{j}]\n")

        # exogenous variables
        f.write(f"\n{indent}# Evaluate numerical values for exogenous variables\n")
        for j, name in enumerate(exo_names):
            f.write(f"{indent}{name} = exo_vars[{j}]\n")

        # parameters
        f.write(f"\n{indent}# Evaluate numerical values for parameters from params\n")
        for j, name in enumerate(param_names):
            f.write(f"{indent}{name} = params[{j}]\n")

        # steady-state variables for dynamic model
        if not is_static:
            f.write(f"\n{indent}# Evaluate numerical values for steady-state variables\\n")
            for j, name in enumerate(endo_names):
                f.write(f"{indent}{name}_ss = steady_state[{j}]\n")

        # initialize output array
        n_rows, n_cols = Output.shape
        f.write(f"\n{indent}# Initialize \n")
        f.write(f"{indent}{name_of_output} = [[0.0]*{n_cols} for _ in range({n_rows})]\n")

        # write non-zero entries only
        f.write(f"\n{indent}# Evaluate non-zero entries\n")
        Output_simplified = Output.applyfunc(sp.simplify)
        for i in range(n_rows):
            for j in range(n_cols):
                expr = Output_simplified[i, j]
                if expr != 0:
                    f.write(f"{indent}{name_of_output}[{i}][{j}] = {sp.python(expr)}\n")

        # return statement
        f.write(f"\n{indent}return {name_of_output}\n")
