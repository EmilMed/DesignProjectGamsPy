from gamspy import (Container, Set, Parameter, Variable,
                    Equation, Model, Sense, Problem, Sum)
import pandas as pd

import numpy as np

# Define model container
m = Container()

# ===============================================================================#
#                         || Base Information ||
# ===============================================================================#

i = Set(
    container=m,
    name='comps',
    records=[
        'H',
        'H2O',
        'O',
        'N',
        'Ar',
        'CO2',
        'NH3'
    ],
    description="Involved chemical components"
)

# No. of streams in process
N_streams = 11

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

# ===============================================================================#
#                               || Crossover 1 ||
# ===============================================================================#

yHfeed = Parameter(
    container=m,
    name='yHfeed',
    domain=[i],
    records=[
        ('H', 0.999070),
        ('H2O', 0.000860),
        ('O', 0.000070),
        ('N', 0.0),
        ('Ar', 0.0),
        ('CO2', 0.0),
        ('NH3', 0.0)
    ],
    description="Molar fraction of components in fresh H2 feed"
)

yNfeed = Parameter(
    container=m,
    name='yNfeed',
    domain=[i],
    records=[
        ('H', 0.0),
        ('H2O', 0.000034),
        ('O', 0.001940),
        ('N', 0.99002),
        ('Ar', 0.007970),
        ('CO2', 0.000036),
        ('NH3', 0.0)
    ],
    description="Molar fraction of components in fresh N2 feed"
)

Crossover1MB = Equation(
    container=m,
    name="Crossover1MB",
    domain=[i],  # over components
    description="Mass balance over crossover 1"
)
Crossover1MB[i] = F[1, i] + F[2, i] == F[3, i]


Hfeed_relation = Equation(
    container=m,
    name="HFeed_relation",
    domain=[i],  # over components
    description="Relating flow rates of fed components in H2"
)
Hfeed_relation[i] = F[1, 'H'] == yHfeed[i] / yHfeed['H'] * F[1, i]


Nfeed_relation = Equation(
    container=m,
    name="NFeed_relation",
    domain=[i],  # over components
    description="Relating flow rates of fed components in N2"
)
Nfeed_relation[i] = F[2, 'N'] == yNfeed[i] / yNfeed['N'] * F[2, i]

# ===============================================================================#
#                               || Crossover 2 ||
# ===============================================================================#

Crossover2MB = Equation(
    container=m,
    name="Crossover2MB",
    domain=[i],  # over components
    description="Mass balance over crossover 2"
)
Crossover2MB[i] = F[3, i] + F[9, i] == F[4, i]


# ===============================================================================#
#                               || Reactor  ||
# ===============================================================================#

v_1 = Parameter(
    container=m,
    name='v_1',
    domain=[i],
    records=[
        ('H', -3),
        ('H2O', 0),
        ('O', 0),
        ('N', -1),
        ('Ar', 0),
        ('CO2', 0),
        ('NH3', 2)
    ],
    description="Stochiometric coefficients for z\
                 reactants in reaction 1 (Haber-Bosch)"
)

sp_conv_hb = Parameter(
    container=m,
    name='sp_conv_hb',
    records=0.12,  # 12% conversion in Haber-Bosch reactor
    description="Specified conversion in Haber-Bosch reactor"
)

Limit_reac_HB = "N"

X_HB = Variable(
    container=m,
    name="X_HB",
    type="positive",
    description="Extent of reaction in Haber-Bosch reactor"
)

Extent_Definition_HB = Equation(
    container=m,
    name="Extent_Definition_HB",
    description="Extent definition in Haber-Bosch reactor"
)

Extent_Definition_HB[...] = X_HB == sp_conv_hb * \
                            F[4, Limit_reac_HB] / (-1 * v_1[Limit_reac_HB])

HB_Reactor_MB = Equation(
    container=m,
    name="HB_Reactor_MB",
    domain=[i],  # over components
    description="Mass balance over Haber-Bosch reactor"
)

HB_Reactor_MB[i] = F[5, i] == F[4, i] + v_1[i] * X_HB


# ===============================================================================#
#                               || Flash  ||
# ===============================================================================#

# Component Mass Balance (Stream 5 -> 6 + 7)
OverallFlashMB = Equation(
    container=m,
    name="OverallFlashMB",
    domain=[i], 
    description="Mass balance over Flash Separator"
)
OverallFlashMB[i] = F[5, i] == F[6, i] + F[7, i]

# Purity bounds for Ammonia!
NH3_Purity_Min = Parameter(container=m,
                           name='Purity_Min',
                           records=0.990,
                           description="Minimum allowed liquid NH3 purity"
                           )

NH3_Purity_Max = Parameter(container=m,
                           name='Purity_Max',
                           records=0.996,
                           description="Maximum allowed liquid NH3 purity")

# Purity Inequality Equations

# NH3_Purity_Lower_Bound = Equation(container=m,
#                                   name='Purity_Lower_Bound',
#                                   description="Liquid NH3 purity must be >="
#                                   )
# NH3_Purity_Lower_Bound[...] = F[6, 'NH3'] / Sum(i, F[6, i]) >= NH3_Purity_Min

# NH3_Purity_Upper_Bound = Equation(container=m,
#                                   name='Purity_Upper_Bound',
#                                   description="Liquid NH3 purity must be <="
#                                   )
# NH3_Purity_Upper_Bound[...] = F[6, 'NH3'] / Sum(i, F[6, i]) <= NH3_Purity_Max


# NH3 Recovery Definition
NH3_recovery = Parameter(
    container=m,
    name="NH3_recovery",
    records=0.98,
    description="NH3 recovery in flash unit"
)

NH3_Recov_Definition = Equation(
    container=m,
    name="NH3_Recov_Definition",
    description="Mass balance over Flash Separator"
)
NH3_Recov_Definition[...] = F[6, 'NH3'] == NH3_recovery*F[5, 'NH3']


i_gases_only = Set(container=m, 
                   name='i_gases_only', 
                   domain=[i], 
                   records=['H', 'O', 'N', 'Ar', 'CO2'],
                   description="Permanent gases forced to vapor phase")

Zero_Gas_Flow = Equation(
    container=m, 
    name='Zero_Gas_Flow',
    domain=[i_gases_only],
    description="Permanent gases forced to vapor phase (F[6,i] == 0)"
)
Zero_Gas_Flow[i_gases_only] = F[6, i_gases_only] == 0


H2O_recovery = Variable(
    container=m,
    name="H2O_recovery",
    type="positive",
    description="H2O recovery in flash unit"
)

H2O_Recov_Definition = Equation(
    container=m,
    name="H2O_Recov_Definition",
    description="H2O definition in flash unit"
)
H2O_Recov_Definition[...] = F[6, 'H2O'] == H2O_recovery*F[5, 'H2O']



# ===============================================================================#
#                               || Split 1  ||
# ===============================================================================#

Split1MB = Equation(
    container=m,
    name="Split1MB",
    domain=[i],  # over components
    description="Mass balance over split 1"
)
Split1MB[i] = F[7, i] == F[8, i] + F[9, i]

HB_PurgeFraction = Parameter(container=m,
                             name='HB_PurgeFraction',
                             records=0.05,
                             description="Purge fraction from HB recycle"
                             )


HB_PurgeFrac_def = Equation(
                            container=m,
                            domain=[i],
                            name="HB_PurgeFrac_def",
                            description="Definition of purge \
                                         fraction from HB recycle"
                            )
HB_PurgeFrac_def[i] = F[8, i] == HB_PurgeFraction * F[7, i]


# ===============================================================================#
#                           || Maritime Transport ||
# ===============================================================================#

MT_MB = Equation(
    container=m,
    name="MT_MB",
    domain=[i],  # over components
    description="Mass balance over MTransport"
)
MT_MB[i] = F[6, i] == F[10, i] + F[11, i]

split_ure = Parameter(
    container=m,
    name='split_ure',
    records=0.6,  # 60% of ammonia to urea process
    description="Split fraction of ammonia from transport to urea process"
)

MT_Split_def = Equation(
                            container=m,
                            domain=[i],
                            name="MT_Split_def",
                            description="MT Split to urea"
                            )
MT_Split_def[i] = F[10, i] == split_ure * F[6, i]


# ===============================================================================#
#                           || PRODUCTION TARGET ||
# ===============================================================================#

NH3_required_py = 46000  # ton/yr
NH3_MR = 17.031  # g/mol
Years_operated = 15  # yr
hours_per_year = 8000  # hr/yr

NH3_mass_production = NH3_required_py * 1e3 / (Years_operated*hours_per_year)  # kg/yr
NH3_molar_production = NH3_mass_production / NH3_MR  # kmol/yr


# Parameter to hold the fixed production rate
NH3_Prod_Target = Parameter(container=m,
                            name='NH3_Prod_Target',
                            records=NH3_molar_production,
                            description="Target NH3 production rate [kmol/hr] in Stream 6")

# Equation to fix the molar flow of NH3 in stream 6
Prod_Target_Constraint = Equation(
    container=m,
    name='Prod_Target_Constraint',
    description="Fixes NH3 output flowrate in Stream 6"
)

# F[6, 'NH3'] must equal the target amount.
Prod_Target_Constraint[...] = F[6, 'NH3'] == NH3_Prod_Target


# ===============================================================================#
#                           || INITIALIZATION FIX ||
# ===============================================================================#

# 1. Initialize all flows F to a small non-zero value for safety
F.l[j, i] = 1e-6  # Small non-zero initial guess for all flows

# 2. Provide a better initial guess for the key variables (Stream 6 and Recoveries)
# You need a non-zero value for the total flow in Stream 6.
# Use the fixed production target (22.5085 kmol/hr, assuming calculation is correct)

# # Set the flow of NH3 in stream 6 to the target value (best guess)
F.lo[6, 'NH3'] = NH3_molar_production 

# # Set the flow of H2O in stream 6 (the main impurity) to a small, non-zero value (e.g., 0.5% of the NH3 flow)
F.lo[6, 'H2O'] = F.l[6, 'NH3'] * 0.005

# # Set initial recovery variables to a non-zero, plausible value (e.g., 0.5)
# NH3_recovery.lo[...] = 0.1
# NH3_recovery.up[...] = 1.0

H2O_recovery.lo[...] = 1e-6
H2O_recovery.up[...] = 1.0


# ===============================================================================#
#                            || MODEL SETUP AND SOLVE ||
# ===============================================================================#

z = Variable(
    container=m,
    name="objectiveZ",
    description="Objective Function Variable")

ObjFunc = Equation(
    container=m,
    name="ObjFunc",
    description="Objective Function Definition"
)

ObjFunc[...] = z == F[1, 'H'] + F[2, 'N']    # Minimize fresh H2 feed


# Define the Model
HB_Process_Model = Model(
    container=m,
    name='Haber_Bosch_Process',
    
    # 1. Set the objective to the H2 Fresh Feed flow
    objective=z,
    sense=Sense.MIN,
    
    # 2. List all required equations (Must include the newly added MT_Split_rem)
    equations=m.getEquations(),
    
    # Use NLP because of the non-linear purity constraints (F[6,NH3]/Sum(F[6,i]))
    problem=Problem.NLP 
)

# Solve the Model
print(HB_Process_Model.solve())
# Display results
print(F.records)
