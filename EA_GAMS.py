from gamspy import (Container, Set, Parameter, Variable,
                    Equation, Model, Sense, Problem, Sum)
import pandas as pd
import numpy as np

# Define model container
m = Container()

# ===============================================================================#
#                         || Base Information ||
# ===============================================================================#

# Involved chemical components 
i = Set(
    container=m,
    name='comps',
    records=[
        'NH3',  'H2O',  'EO',   'MEA',  'DEA',  'TEA',  'Heavy', 
        'H',    'O',    'N',    'CO2',  'Ar'
    ],
    description="Involved chemical components including inert impurities"
)

# No. of streams in process
N_streams_start = 23
N_streams_end = 39
j = Set(
    container=m,
    name='streams',
    records=list(range(N_streams_start, N_streams_end + 1)),
    description="Flow streams in process"
)

# Molar flowrate of component i in stream j
F = Variable(
    container=m,
    name="F",
    domain=[j, i],
    type="positive",
    description="Molar flowrate of component i in stream j [kmol/hr]"
)

# Extent of reaction variable
X = Variable(
    container=m,
    name="X",
    type="positive",
    description="Extent of reaction for the PFR Reactor [kmol/hr]"
)

# Split Fraction for Splitter (34 -> 32 + 35)
sf = Variable(
    container=m, 
    name='sf', 
    type='positive', 
    description="Split Fraction to Purge (Stream 35)"
)

# ===============================================================================#
#                         || STOICHIOMETRY & FEEDS ||
# ===============================================================================#

sp_conv = Parameter(container=m, name='sp_conv', records=0.99, description="Specified EO conversion")
sel_MEA = Parameter(container=m, name='sel_MEA', records=0.75, description="Selectivity to MEA")
sel_DEA = Parameter(container=m, name='sel_DEA', records=0.21, description="Selectivity to DEA")
sel_TEA = Parameter(container=m, name='sel_TEA', records=0.04, description="Selectivity to TEA")

# --- Impure Ammonia Feed (Stream 23) ---
yNH3_feed = Parameter(
    container=m,
    name='yNH3_feed',
    domain=[i],
    records=[
        ('H', 0.0), ('H2O', 0.001784), ('O', 0.000330), ('N', 0.001881), 
        ('CO2', 0.000005), ('NH3', 0.996000), ('Ar', 0.0), ('EO', 0.0),
        ('MEA', 0.0), ('DEA', 0.0), ('TEA', 0.0), ('Heavy', 0.0)
    ]
)

TotalFlow23 = Variable(container=m, name="TotalFlow23", type="positive")
TotalFlow23_Def = Equation(container=m, name="TotalFlow23_Def")
TotalFlow23_Def[...] = TotalFlow23 == Sum(i, F[23, i])

NH3Feed_relation = Equation(container=m, name="NH3Feed_relation", domain=[i])
NH3Feed_relation[i] = F[23, i] == TotalFlow23 * yNH3_feed[i]

# ===============================================================================#
#                         || MIXER (23 + 24 -> 25) ||
# ===============================================================================#

Mixer25MB = Equation(
    container=m,
    name="Mixer25MB",
    domain=[i]
)

Mixer25MB[i] = F[23, i] + F[24, i] == F[25, i]

WaterRatio = Equation(
    container=m,
    name="WaterRatio"
)

WaterRatio[...] = F[25, 'H2O'] == 0.25 * F[25, 'NH3']

# ===============================================================================#
#                           || MIXER (25 + 26 -> 27) ||
# ===============================================================================#

Mixer27MB = Equation(
    container=m,
    name="Mixer27MB",
    domain=[i]
)

Mixer27MB[i] = F[25, i] + F[26, i] == F[27, i]

# ===============================================================================#
#                        || PFR Reactor (27 + 28 -> 29) ||
# ===============================================================================#

ExtentDef = Equation(
    container=m,
    name="ExtentDef"
)

ExtentDef[...] = X == sp_conv * (F[28, 'EO'] + F[27, 'EO'])

ReactorMB = Equation(
    container=m,
    name="ReactorMB"
)

ReactorMB['NH3'] = F[29, 'NH3'] == F[27, 'NH3'] - X
ReactorMB['EO'] = F[29, 'EO'] == F[27, 'EO'] + F[28, 'EO'] - X
ReactorMB['MEA'] = F[29, 'MEA'] == X * sel_MEA
ReactorMB['DEA'] = F[29, 'DEA'] == X * sel_DEA
ReactorMB['TEA'] = F[29, 'TEA'] == X * sel_TEA

InertComps = Set(
    container=m,
    name='InertComps',
    records=['H2O', 'Heavy', 'H', 'O', 'N', 'CO2', 'Ar']
)

ReactorMB[InertComps] = F[29, InertComps] == F[27, InertComps] + F[28, InertComps]

# ===============================================================================#
#                       || Ammonia Stripper (29 -> 30 + 31) ||
# ===============================================================================#

StripperMB = Equation(
    container=m,
    name="StripperMB",
    domain=[i]
)

StripperMB[i] = F[29, i] == F[30, i] + F[31, i]

NH3Recovery = Equation(
    container=m,
    name="NH3Recovery"
)

NH3Recovery[...] = F[30, 'NH3'] == 0.99 * F[29, 'NH3']

NonVolatileSet = Set(
    container=m,
    name='NonVolatileSet',
    records=['H2O', 'MEA', 'DEA', 'TEA', 'Heavy', 'H', 'O', 'N', 'CO2', 'Ar']
)

NonVolatileSplit = Equation(
    container=m,
    name="NonVolatileSplit", 
    domain=[NonVolatileSet]
)

NonVolatileSplit[NonVolatileSet] = F[30, NonVolatileSet] == 0

# ===============================================================================#
#                          || Mixer (30 + 32 -> 26) ||
# ===============================================================================#

Mixer26MB = Equation(
    container=m,
    name="Mixer26MB",
    domain=[i]
)

Mixer26MB[i] = F[30, i] + F[32, i] == F[26, i]

# ===============================================================================#
#                     || Dehydration Unit (31 -> 33 + 34) ||
# ===============================================================================#

DehyMB = Equation(
    container=m,
    name="DehyMB",
    domain=[i]
)

DehyMB[i] = F[31, i] == F[33, i] + F[34, i]

H2ORecovery = Equation(
    container=m,
    name="H2ORecovery"
)

H2ORecovery[...] = F[34, 'H2O'] == 0.99 * F[31, 'H2O']

DehyNonVolatileSet = Set(
    container=m,
    name='DehyNonVolatileSet',
    records=['MEA', 'DEA', 'TEA', 'Heavy', 'H', 'O', 'N', 'CO2', 'Ar']
)

NonVolatileDehySplit = Equation(
    container=m,
    name="NonVolatileDehySplit",
    domain=[DehyNonVolatileSet]
)

NonVolatileDehySplit[DehyNonVolatileSet] = F[34, DehyNonVolatileSet] == 0

# ===============================================================================#
#                      || Splitter (34 -> 32 + 35) ||
# ===============================================================================#

SplitterMB = Equation(
    container=m,
    name="SplitterMB",
    domain=[i]
)

SplitterMB[i] = F[34, i] == F[32, i] + F[35, i]

SplitterFlow35 = Equation(
    container=m,
    name="SplitterFlow35",
    domain=[i]
)

SplitterFlow35[i] = F[35, i] == sf * F[34, i]

# ===============================================================================#
#                  || EAs Separation Unit (33 -> 36, 37, 38, 39) ||
# ===============================================================================#

SeparationMB = Equation(
    container=m,
    name="SeparationMB",
    domain=[i]
)

SeparationMB[i] = F[33, i] == F[36, i] + F[37, i] + F[38, i] + F[39, i]

# Purity and Recovery Constraints (DOF Table)
MEA_Purity = Equation(
    container=m,
    name="MEA_Purity"
)

MEA_Purity[...] = F[36, 'MEA'] == 0.99 * Sum(i, F[36, i])

MEA_Recovery = Equation(
    container=m,
    name="MEA_Recovery"
)

MEA_Recovery[...] = F[36, 'MEA'] == 0.997 * F[33, 'MEA']

DEA_Purity = Equation(
    container=m,
    name="DEA_Purity"
)

DEA_Purity[...] = F[37, 'DEA'] == 0.99 * Sum(i, F[37, i])

DEA_Recovery = Equation(
    container=m,
    name="DEA_Recovery"
)

DEA_Recovery[...] = F[37, 'DEA'] == 0.997 * F[33, 'DEA']

TEA_Purity = Equation(
    container=m,
    name="TEA_Purity"
)

TEA_Purity[...] = F[38, 'TEA'] == 0.99 * Sum(i, F[38, i])

TEA_Recovery = Equation(
    container=m,
    name="TEA_Recovery"
)

TEA_Recovery[...] = F[38, 'TEA'] == 0.997 * F[33, 'TEA']


# --- PRODUCTION TARGET ---
MEA_Prod_Target = Parameter(container=m, name='MEA_Prod_Target', records=100.0)
Prod_Target_Constraint = Equation(container=m, name='Prod_Target_Constraint')
Prod_Target_Constraint[...] = F[36, 'MEA'] == MEA_Prod_Target

# ===============================================================================#
#                           || INITIALIZATION AND BOUNDS ||
# ===============================================================================#

F.l[j, i] = 1e-6 # Initialize all flows
X.l = 100.0
sf.l = 0.05

# Define lists of components that MUST have zero flow in specific streams
# Components that must be ZERO in Stream 24 (Pure H2O feed)
zero_comp_24 = ['NH3', 'EO', 'MEA', 'DEA', 'TEA', 'Heavy', 'H', 'O', 'N', 'CO2', 'Ar'] 
# Components that must be ZERO in Stream 28 (EO/H2O feed)
zero_comp_28 = ['NH3', 'MEA', 'DEA', 'TEA', 'Heavy', 'H', 'O', 'N', 'CO2', 'Ar'] 
# Components that must be ZERO in Stream 37 (DEA product + allowed trace TEA, so zero everything else)
zero_comp_37 = ['NH3', 'H2O', 'EO', 'MEA', 'Heavy', 'H', 'O', 'N', 'CO2', 'Ar'] 
# Components that must be ZERO in Stream 38 (TEA product only)
zero_comp_38 = ['NH3', 'H2O', 'EO', 'MEA', 'DEA', 'Heavy', 'H', 'O', 'N', 'CO2', 'Ar'] 
# Components that must be ZERO in Stream 39 (Heavy Ends only)
zero_comp_39 = ['NH3', 'H2O', 'EO', 'MEA', 'DEA', 'TEA', 'H', 'O', 'N', 'CO2', 'Ar'] 
# Components that must be ZERO in Stream 36 (MEA product + impurities, but DEA/TEA/Heavy should be minimal)
# The zeroing of DEA/TEA/Heavy in 36 is indirectly handled by the high recoveries/purities of 37, 38, 39, 
# so we will not enforce explicit F.up=0 here to avoid over-constraining the model.

# Enforce Purity/Zero Flow via F.up bounds (Robust Method)
def enforce_zero_flow(stream_index, component_list):
    for comp in component_list:
        F.up[stream_index, comp] = 0.0
        F.l[stream_index, comp] = 0.0

enforce_zero_flow(24, zero_comp_24)
enforce_zero_flow(28, zero_comp_28)
enforce_zero_flow(37, zero_comp_37)
enforce_zero_flow(38, zero_comp_38)
enforce_zero_flow(39, zero_comp_39)


# ===============================================================================#
#                            || MODEL SETUP AND SOLVE ||
# ===============================================================================#

z = Variable(container=m, name="objectiveZ")

ObjFunc = Equation(container=m, name="ObjFunc")
# Minimize fresh EO feed (Stream 28)
ObjFunc[...] = z == F[28, 'EO']


# Define the Model
EA_Process_Model = Model(
    container=m,
    name='Ethanolamines_Process',
    objective=z,
    sense=Sense.MIN,
    equations=m.getEquations(),
    problem=Problem.NLP
)

