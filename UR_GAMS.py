from gamspy import (Container, Set, Parameter, Variable,
                    Equation, Model, Sense, Problem, Sum)
import pandas as pd


def null_comp_check(param: Parameter, i: Set):
    """
    Checks if all components in set 'i' are present in the parameter's
    records, and adds them with a value of 0.0 if missing.
    """
    for component in i.records.iloc[:, 0]:
        if component not in [record[0] for record in param.records]:
            param.records.loc[len(param.records)] = [component, 0.0]
    return param


def fix_values(var: Variable, val: float):
    """Fixes the lower, upper, and level bounds \
        of a Variable to a specific value."""
    var.up[...] = val
    var.l[...] = val
    var.lo[...] = val
    return var


# Define model container
m = Container()

# ===============================================================================#
#                          || Base Information ||
# ===============================================================================#

i = Set(
    container=m,
    name='comps',
    records=[
        'H',
        'H2O',
        'O',
        'N',
        'CO2',
        'NH3',
        'Ar',
    ],
    description="Involved chemical components"
)

# No. of streams in process
N_streams = -

j = Set(
    container=m,
    name='streams',
    records=list(range(1, N_streams + 1)),
    description="Flow streams in process"
)


F = Variable(
    container=m,
    name="F",
    domain=[j, i],  # per comp & stream
    type="positive",
    description="Molar flowrate of component i in stream j"
)