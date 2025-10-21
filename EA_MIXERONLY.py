from gamspy import (Container, Set, Parameter, Variable,
                    Equation, Model, Sense, Problem, Sum)
import pandas as pd


def null_comp_check(param: Parameter, i: Set):
    """
    Checks if all components in set 'i' are present in the parameter's
    records, and adds them with a value of 0.0 if missing.
    """
    # This function is not used in this specific block but is kept for completeness
    for component in i.records.iloc[:, 0]:
        if component not in [record[0] for record in param.records]:
            param.records.loc[len(param.records)] = [component, 0.0]
    return param


def fix_values(var: Variable, val: float):
    """Fixes the lower, upper, and level bounds 
        of a Variable to a specific value."""
    var.up[...] = val
    var.l[...] = val
    var.lo[...] = val
    var.fx[...] = val
    return var


# Define model container
m = Container()

# ===============================================================================#
#                              || Base Information ||
# ===============================================================================#

i = Set(
    container=m,
    name='comps',
    records=[
        'NH3',
        'H2O',
        'EO',
        'MEA',
        'DEA',
        'TEA',
    ],
    description="Involved chemical components EA process"
)

# No. of streams in process
N_streams_start = 23
N_streams_end = 25

j = Set(
    container=m,
    name='streams',
    records=list(range(N_streams_start, N_streams_end + 1)),
    description="Flow streams in process"
)


F = Variable(
    container=m,
    name="F",
    domain=[j, i],  # per comp & stream
    type="positive",
    description="Molar flowrate of component i in stream j"
)

# --- Initialization and Fixed Inputs ---
# F.lo[j, i] = 1e-6
# F.l[j, i] = 0.5  # Initial guess for solver stability

# ===============================================================================#
#                              || Mixer 1 (23+24->25)||
# ===============================================================================#

Mixer1mb = Equation(
    container=m,
    domain=i,
)
Mixer1mb[i] = F[25, i] == F[23, i] + F[24, i]    # 6 equations

notwater = Set(
    container=m,
    domain = i,
    name='notwater',
    records=[
        'NH3',
        'EO',
        'MEA',
        'DEA',
        'TEA'
    ],
    description="Involved chemical components EA process excl water"
)

S24_Comp = Equation(
    container=m,
    domain=notwater,
)
S24_Comp[notwater] = F[24, notwater] == 0.0  # 5 equations

NH3_input = 135.015  # kmol/hr
fix_values(F[23, 'NH3'], NH3_input)
fix_values(F[23, 'H2O'], 0.0)
fix_values(F[23, 'EO'], 0.0)
fix_values(F[23, 'MEA'], 0.0)
fix_values(F[23, 'TEA'], 0.0)
fix_values(F[23, 'DEA'], 0.0)  # 6 equations


WaterComp = Equation(
    container=m
)
WaterComp[...] = F[25, 'H2O']*3 == F[25, 'NH3']  # 1 equation


#===============================================================================#
#                         || OBJECTIVE FUNCTION ||
# ===============================================================================#

F24H2O_obj = Variable(
    container=m,
    name="F24H2O_obj",
    description="Objective: Minimize H2O Flow Rate in Stream 24"
)

ObjDef = Equation(
    container=m,
    name="ObjDef",
    description="Defines objective variable as F[24, H2O]"
)
ObjDef[...] = F24H2O_obj == F[24, 'H2O']

# ===============================================================================#
#                              || MODEL SETUP AND SOLVE ||
# ===============================================================================#

Mixer_Model = Model(
    container=m,
    name='Mixer_Block_Model',
    # Equations are now 6 (MB) + 5 (S24_Comp) + 1 (WaterComp) + 1 (ObjDef) = 13 total equations
    # The objective adds a new variable (F24H2O_obj) and an equation (ObjDef)
    equations=[Mixer1mb, S24_Comp, WaterComp, ObjDef], 
    objective=F24H2O_obj,
    sense=Sense.MIN,
    problem=Problem.NLP
)

print(Mixer_Model.solve())