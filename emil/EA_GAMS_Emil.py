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
    # var.up[...] = val
    # var.l[...] = val
    # var.lo[...] = val
    var.fx[...] = val
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
N_streams_end = 39

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


# ===============================================================================#
#                              || Mixer 1 (23+24->25)|| W
# ===============================================================================#

Mixer1mb = Equation(
    container=m,
    domain=i,
)
Mixer1mb[i] = F[23, i] + F[24, i] == F[25, i]

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
S24_Comp[notwater] = F[24, notwater] == 0.0

WaterComp = Equation(
    container=m,
)
WaterComp[...] = F[25, 'H2O']*3 == F[25, 'NH3']

# ===============================================================================#
#                              || Crossover 1 (25+26->27)|| W
# ===============================================================================#

Crossover1mb = Equation(
    container=m,
    domain=i,
)
Crossover1mb[i] = F[25, i] + F[26, i] == F[27, i]

# ===============================================================================#
#                              ||   PFR 1  (27+28->29)||
# ===============================================================================#
reactions = Set(
    container=m,
    name='reactions',
    records=[1, 2, 3],
    description="Reactions in PFR"
)

nu = Parameter(
    container=m,
    domain=[i, reactions],
    name='nu',
    records=[
        ('NH3', 1, -1.0),
        ('NH3', 2, 0),
        ('NH3', 3, 0),

        ('H2O', 1, 0.0),
        ('H2O', 2, 0.0),
        ('H2O', 3, 0.0),

        ('EO', 1, -1.0),
        ('EO', 2, -1.0),
        ('EO', 3, -1.0),

        ('MEA', 1, 1),
        ('MEA', 2, -1.0),
        ('MEA', 3, 0),

        ('DEA', 1, 0),
        ('DEA', 2, 1.0),
        ('DEA', 3, -1.0),

        ('TEA', 1, 0),
        ('TEA', 2, 0),
        ('TEA', 3, 1.0),
    ],
    description="Stoichiometric coefficients for the PFR"
)

SPC_EO = Parameter(
    container=m,
    name='SPC_EO',
    records=0.99
)

S_MEA = Parameter(
    container=m,
    name='S_MEA',
    records=0.75
)

S_DEA = Parameter(
    container=m,
    name='S_DEA',
    records=0.21
)

# blahhh
S_TEA = Parameter(
    container=m,
    name='S_TEA',
    records=0.04
)

# MAYBE
extent = Variable(
    container=m,
    name="extent_of_reaction",
    domain=reactions,
    type="positive",
    description="Extent of reaction for each reaction in PFR"
)

pfr1mb = Equation(
    container=m,
    domain=i,
)
pfr1mb[i] = F[27, i] + F[28, i] + Sum(reactions, nu[i, reactions] * extent[reactions]) == F[29, i]

SPCdef = Equation(
    container=m,
)
SPCdef[...] = Sum(reactions, extent[reactions]) == SPC_EO * (F[27, 'EO'] + F[28, 'EO'])

S_MEAdef = Equation(
    container=m,
)
S_MEAdef[...] = (extent[1]-extent[2]) == (extent[1]+extent[2]+extent[3]) * S_MEA

S_DEAdef = Equation(
    container=m,
)
S_DEAdef[...] = 2*(extent[2]-extent[3]) == (extent[1]+extent[2]+extent[3]) * S_DEA

SumEAs = Equation(
    container = m,
)
SumEAs[...] = SPC_EO * (F[27, 'EO'] + F[28, 'EO']) == F[29, 'MEA'] + F[29, 'DEA'] + F[29, 'TEA']

NH3EOratio = Equation(
    container=m,
)
NH3EOratio[...] = (F[27, 'EO'] + F[28, 'EO']) * 10 == F[27, 'NH3']

notwaterANDEO = Set(
    container=m,
    domain=i,
    name='notwaterANDEO',
    records=[
        'NH3',
        'MEA',
        'DEA',
        'TEA'
    ],
    description="Involved chemical components EA process excl water and eo"
) 

S28_Comp = Equation(
    container=m,
    domain=notwaterANDEO,
)
S28_Comp[notwaterANDEO] = F[28, notwaterANDEO] == 0.0

S28FeedRatio = Equation(
    container=m,
)

S28FeedRatio[...] = F[28, 'EO'] == 99 * (F[28, 'H2O'])

# ===============================================================================#
#                              ||   Ammonia Stripper (29->30+31)||
# ===============================================================================#

Strippermb = Equation(
    container=m,
    domain=i,
)

Strippermb[i] = F[29, i] == F[30, i] + F[31, i]

ammoniaRecovery = Equation(
    container=m,
)

ammoniaRecovery[...] = F[30, 'NH3'] == 0.99 * F[29, 'NH3']

notAmmonia = Set(
    container=m,
    domain=i,
    name='notAmmonia',
    records=[
        'MEA',
        'DEA',
        'TEA',
        'H2O',
        'EO'
    ],
    description="Involved chemical components EA process excl ammonia"
)

otherRecovery = Equation(
    container=m,
    domain=notAmmonia,
)

otherRecovery[notAmmonia] = F[31, notAmmonia] == F[29, notAmmonia]

# ===============================================================================#
#                              ||   Crossover 2  (30+32->26) ||
# ===============================================================================# 

Crossover2mb = Equation(
    container=m,
    domain=i,
)

Crossover2mb[i] = F[30, i] + F[32, i] == F[26, i]


# ===============================================================================#
#                              ||   Dehydration Unit (31->33+34) ||
# ===============================================================================# 

Dehydrationmb = Equation(
    container=m,
    domain=i,
)

Dehydrationmb[i] = F[31, i] == F[33, i] + F[34, i]

waterRemoval = Equation(
    container=m,
)
waterRemoval[...] = F[33, 'H2O'] == 0.99 * F[31, 'H2O']

notNH3EOH2O = Set(
    container=m,
    domain=i,
    name='notNH3EOH2O',
    records=[
        'MEA',
        'DEA',
        'TEA',
    ],
    description="Involved chemical components EA process excl ammonia"
)

otherRecDehyd = Equation(
    container=m,
    domain=notNH3EOH2O,
)
otherRecDehyd[notNH3EOH2O] = F[31, notNH3EOH2O] == F[34, notNH3EOH2O]

DehydEO = Equation(
    container=m,
)
DehydEO[...] = F[31, 'EO'] == F[33, 'EO']

DehydNH3 = Equation(
    container=m,
)
DehydNH3[...] = F[31, 'NH3'] == F[33, 'NH3']

# ===============================================================================#
#                              ||   Splitter (33->32+35) ||
# ===============================================================================# 

sf = Parameter(
    container=m,
    name='sf',
    records=0.99,
)

SplitterMB = Equation(
    container=m,
    domain=i,
)
SplitterMB[i] = F[32, i] + F[35, i] == F[33, i]

Splitratio = Equation(
    container=m,
    domain=i,
)
Splitratio[i] = sf * F[33, i] == F[35, i]

# ===============================================================================#
#                              ||   Distillation Column (34->36+37+38+39) ||
# ===============================================================================# 

RecoveryDMEA = Parameter(
    container=m,
    name='RecoveryDMEA',
    records = 0.999999
)

RecoveryDDEA = Parameter(
    container=m,
    name='RecoveryDDEA',
    records = 0.999999
)

RecoveryDTEA = Parameter(
    container=m,
    name='RecoveryDTEA',
    records = 0.999999
)

DistillationMB = Equation(
    container=m,
    domain=i,
)

DistillationMB[i] = F[34, i] == F[36, i] + F[37, i] + F[38, i] + F[39, i]

# NotH2OEA = Set(
#     container=m,
#     domain=i,
#     name='NotH2OEA',
#     records=[
#         'NH3',
#         'EO'
#     ],
#     description="Involved chemical components excl EAs and water"
# )

# RECOVERIES
RecoveryH2ODef = Equation(
    container=m,
)
RecoveryH2ODef[...] = F[36, 'H2O'] == 1.0*F[34, 'H2O']

RecoveryMEADef = Equation(
    container=m,
)
RecoveryMEADef[...] = F[36, 'MEA'] == RecoveryDMEA*F[34, 'MEA']

RecoveryDEADef = Equation(
    container=m,
)
RecoveryDEADef[...] = F[37, 'DEA'] == RecoveryDDEA*F[34, 'DEA']

RecoveryTEADef = Equation(
    container=m,
)
RecoveryTEADef[...] = F[38, 'TEA'] == RecoveryDTEA*F[34, 'TEA']

# 1 - RECOVERIES
RecoveryMEA2Def = Equation(
    container=m,
)
RecoveryMEA2Def[...] = F[37, 'MEA'] == (1-RecoveryDMEA)*F[34, 'MEA']

RecoveryDEA2Def = Equation(
    container=m,
)
RecoveryDEA2Def[...] = F[38, 'DEA'] == (1-RecoveryDDEA)*F[34, 'DEA']

RecoveryTEA3Def = Equation(
    container=m,
)
RecoveryTEA3Def[...] = F[39, 'TEA'] == (1-RecoveryDTEA)*F[34, 'TEA']   

# # PURITIES
# PurityDMEA = Equation(
#     container=m,
# )
# PurityDMEA[...] = F[36, 'MEA'] == 0.99*Sum(i, F[36, i])

# PurityDDEA = Equation(
#     container=m,
# )

# PurityDDEA[...] = F[37, 'DEA'] == 0.99*Sum(i, F[37, i])

# PurityDTEA = Equation(
#     container=m,
# )
# PurityDTEA[...] = F[38, 'TEA'] == 0.99*Sum(i, F[38, i])

# ===============================================================================#
#                             || INLET REQUIREMENTS ||
# ===============================================================================#
# ===============================================================================#
#                             || INITIALIZATION FIX ||
# ===============================================================================#
# 1. Initialize all flows F to a small non-zero value for safety
F.l[j, i] = 1e-6  # Small non-zero initial guess for all flows

NH3_input = 135.015  # kmol/hr
fix_values(F[23, 'NH3'], NH3_input)
fix_values(F[23, 'H2O'], 0.0)
fix_values(F[23, 'EO'], 0.0)
fix_values(F[23, 'MEA'], 0.0)
fix_values(F[23, 'TEA'], 0.0)
fix_values(F[23, 'DEA'], 0.0)

# Set initial recovery variables
# RecoveryDMEA.lo[...] = 0.8
# RecoveryDMEA.up[...] = 0.9

# RecoveryDDEA.lo[...] = 0.8
# RecoveryDDEA.up[...] = 0.9

# RecoveryDTEA.lo[...] = 0.8
# RecoveryDTEA.up[...] = 0.9

# ===============================================================================#
#                            || MODEL SETUP AND SOLVE ||
# ===============================================================================#

z = Variable(
    container=m,
    name="objectiveZ",
    description="Objective Function Variable"
)

ObjFunc = Equation(
    container=m,
    name="ObjFunc",
    description="Objective Function Definition"
)

ObjFunc[...] = z == Sum(reactions, extent[reactions])

# Define the Model
EA_Prodution_Model = Model(
    container=m,
    name='EA_Prodution_Model',

    # 1. Set the objective to the H2 Fresh Feed flow
    objective=z,
    sense=Sense.MAX,

    # 2. List all required equations
    equations=m.getEquations(),

    # Use NLP because of the non-linear purity
    problem=Problem.NLP
)

print(EA_Prodution_Model.solve())

pd.set_option('display.max_rows', None)     # Display all rows
pd.set_option('display.max_columns', None)  # Display all columns
pd.set_option('display.width', None)        # Allow output to be wider

print("--- Component Flows per Stream ---")
print(F.records)

data = F.records.values

# 1. Create a DataFrame from the raw data
# Assign generic names based on the columns of interest
df = pd.DataFrame(data, columns=['stream_number', 'component', 'flow_level',
                                 'marginal', 'lower', 'upper', 'scale'])

# 2. Ensure the key columns are numeric (critical for summation!)
# Errors='coerce' handles the string values
# like 'inf' by converting them to NaN,
# which are correctly ignored by the .sum() function.

df['stream_number'] = pd.to_numeric(df['stream_number'], errors='coerce')
df['flow_level'] = pd.to_numeric(df['flow_level'], errors='coerce')

# 3. The Formula: Group by stream_number and sum the flow_level
stream_sums = df.groupby('stream_number')['flow_level'].sum().reset_index()

# Rename columns for clarity
stream_sums.rename(columns={'flow_level': 'total_flow_sum'}, inplace=True)

# Print the result
print("--- Total Flow Sums per Stream ---")
print(stream_sums)

# ----------------------------------------------------------------------
# ---           Composition Calculation                            ---

print("\n--- Composition of Streams (Component Fraction) ---")

# 4. Merge the total stream sums back into the original data frame (df).
# The merge works cleanly because stream_sums now only contains
# 'stream_number' and 'total_flow_sum'
composition_df = df.merge(
    stream_sums,
    on='stream_number',
    how='left'
)

# 5. Calculate the fractional composition for each component.
# FIX: The numerator is correctly 'flow_level' (from original df).
# The denominator is 'total_flow_sum' (the new column from the merge).
composition_df['fraction'] = (
    composition_df['flow_level'] / composition_df['total_flow_sum']
)

# Optional: Select the relevant columns
# and round the fraction for cleaner output
final_composition = composition_df[['stream_number', 'component', 'flow_level',
                                    'total_flow_sum', 'fraction']].copy()
final_composition['fraction'] = final_composition['fraction'].round(6)

# Final display
print(final_composition)
