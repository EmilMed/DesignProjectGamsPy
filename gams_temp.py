from gamspy import (Container, Set, Parameter, Variable,
                    Equation, Model, Sum, Sense, Options)

# Define model container
m = Container()

# ===============================================================================#
#                                 || Sets ||
# ===============================================================================#

# Needs to be defined with actual data

i = Set(
    Container=m,
    name='comps',
    records=[
        'H',
        'H2O',
        'O',
        'N',
        'Ar',
        'CO2',
    ],
    description="Involved chemical components"
)

j = Set(Container=m, name='streams', description="Flow streams in process")


# ===============================================================================#
#                              || Parameters ||
# ===============================================================================#

# Needs to be defined with actual data

v1 = Parameter(
    Container=m,
    name='v1',
    domain=i,  # over components
    description="Stoichiometric coefficient of component i in haber reactor"
)

X_ure = Parameter(
    Container=m,
    name='X_Ure',
    description="Split fraction of ammonia from transport to urea process"
)

# ===============================================================================#
#                              || Variables ||
# ===============================================================================#

F = Variable(
    Container=m,
    name="F",
    domain=[i, j],  # per comp & stream
    type="Positive",
    description="Molar flowrate of component i in stream j"
)


X_haber = Parameter(
    Container=m,
    name='X_haber',
    type="Positive",
    description="Split fraction of reactant gas to be recycled from \
                 flash drum to haber be mixed with fresh feed"
)

# ===============================================================================#
#                              || Equations ||
# ===============================================================================#
