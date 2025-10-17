from gamspy import (Container, Set, Parameter, Variable,
                    Equation, Model, Sum, Sense, Options)

# Define model container
m = Container()

# ===============================================================================#
#                            || Base Information ||
# ===============================================================================#


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

# No. of streams in process
N_streams = 23

j = Set(
    Container=m,
    name='streams',
    records=list(range(1, N_streams + 1)),
    description="Flow streams in process"
)

X_ure = Parameter(
    Container=m,
    name='X_Ure',
    records=0.6,  # 60% of ammonia to urea process
    description="Split fraction of ammonia from transport to urea process"
)


F = Variable(
    Container=m,
    name="F",
    domain=[j, i],  # per comp & stream
    type="Positive",
    description="Molar flowrate of component i in stream j"
)

X_haber = Variable(
    Container=m,
    name='X_haber',
    type="Positive",
    description="Split fraction of reactant gas to be recycled from \
                 flash drum to haber be mixed with fresh feed"
)

# ===============================================================================#
#                                  || Crossover 1 ||
# ===============================================================================#

yHfeed = Parameter(
    Container=m,
    name='yHfeed',
    domain=[i],
    records=[
        ('H', 0.999070),
        ('H2O', 0.000860),
        ('O', 0.000070),
        ('N', 0.0),
        ('Ar', 0.0),
        ('CO2', 0.0)
    ],
    description="Molar fraction of components in fresh H2 feed"
)

yNfeed = Parameter(
    Container=m,
    name='yHfeed',
    domain=[i],
    records=[
        ('H', 0.0),
        ('H2O', 0.000034),
        ('O', 0.001940),
        ('N', 0.99002),
        ('Ar', 0.007970),
        ('CO2', 0.000036)
    ],
    description="Molar fraction of components in fresh N2 feed"
)

Crossover1 = Equation(
    Container=m,
    name="Crossover1",
    domain=[i],  # over components
    description="Mass balance over crossover 1"
)[i] = F[1, i] + F[2, i] == F[3, i]

Hfeed_relation = Equation(
    Container=m,
    name="HFeed_relation",
    domain=[i],  # over components
    description="Relating flow rates of fed components in H2"
)[i] = yHfeed['H'] * F[1, 'H'] == yHfeed[i] * F[1, i]

Nfeed_relation = Equation(
    Container=m,
    name="NFeed_relation",
    domain=[i],  # over components
    description="Relating flow rates of fed components in N2"
)[i] = yNfeed['N'] * F[2, 'N'] == yNfeed[i] * F[2, i]


# ===============================================================================#
#                                  || Crossover 2 ||
# ===============================================================================#

Crossover2 = Equation(
    Container=m,
    name="Crossover2",
    domain=[i],  # over components
    description="Mass balance over crossover 2"
)[i] = F[3, i] + F[9, i] == F[4, i]


# ===============================================================================#
#                                  || Reactor ||
# ===============================================================================#

v1 = Parameter(
    Container=m,
    name='v1',
    domain=[i],  # over components
    records=[
        ('H', 3),
        ('H2O', 0),
        ('O', 0),
        ('N', 1),
        ('Ar', 0),
        ('CO2', 0)
    ],
    description="Stoichiometric coefficient of component i in haber reactor"
)