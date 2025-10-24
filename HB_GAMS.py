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
#                           || STOICHIOMETRIC FEED ||
# ===============================================================================#

# BAD COND
 
# # Equation to enforce the H:N ratio for fresh feed (Stream 3)
# StoichFeed = Equation(
#     container=m,
#     name='StoichFeed',
#     description="Fresh feed ratio must be H:N = 3:1"
# )
# # F[1, 'H'] is the H2 component flow. F[2, 'N'] is the N2 component flow.
# StoichFeed[...] = F[1, 'H'] == 3 * F[2, 'N']


# ===============================================================================#
#                              || Crossover 1 ||
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

null_comp_check(yHfeed, i)

yNfeed = Parameter(
    container=m,
    name='yNfeed',
    domain=[i],
    records=[
        # ('H', 0.0),
        ('H2O', 0.000034),
        ('O', 0.001940),
        ('N', 0.99002),
        ('Ar', 0.007970),
        ('CO2', 0.000036),
        # ('NH3', 0.0)
    ],
    description="Molar fraction of components in fresh N2 feed"
)

null_comp_check(yNfeed, i)

TotalFlow1 = Variable(
    container=m,
    name="TotalFlow1",
    type="positive",
    description="Total molar flowrate of fresh H2 feed (Stream 1)"
)
TotalFlow2 = Variable(
    container=m,
    name="TotalFlow2",
    type="positive",
    description="Total molar flowrate of fresh N2 feed (Stream 2)"
)

TotalFlow1_Def = Equation(
    container=m,
    name="TotalFlow1_Def",
    description="Definition of total molar flowrate of fresh H2 feed"
)
TotalFlow1_Def[...] = TotalFlow1 == Sum(i, F[1, i])


TotalFlow2_Def = Equation(
    container=m,
    name="TotalFlow2_Def",
    description="Definition of total molar flowrate of fresh H2 feed"
)
TotalFlow2_Def[...] = TotalFlow2 == Sum(i, F[2, i])

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
Hfeed_relation[i] = TotalFlow1 * yHfeed[i] == F[1, i]


Nfeed_relation = Equation(
    container=m,
    name="NFeed_relation",
    domain=[i],  # over components
    description="Relating flow rates of fed components in N2"
)
Nfeed_relation[i] = TotalFlow2 * yNfeed[i] == F[2, i]

# ===============================================================================#
#                              || Crossover 2 ||
# ===============================================================================#

Crossover2MB = Equation(
    container=m,
    name="Crossover2MB",
    domain=[i],  # over components
    description="Mass balance over crossover 2"
)
Crossover2MB[i] = F[3, i] + F[9, i] == F[4, i]


# ===============================================================================#
#                                || Reactor  ||
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
    description="Stochiometric coefficients for\
                 reactants in reaction 1 (Haber-Bosch)"
)

sp_conv_hb = Parameter(
    container=m,
    name='sp_conv_hb',
    records=0.12,  # 12% conversion in Haber-Bosch reactor
    description="Specified conversion in Haber-Bosch reactor"
)

Limit_reac_HB = "H"  # Limiting reactant for Haber-Bosch reaction

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

Extent_Definition_HB[...] = (
    X_HB == sp_conv_hb * F[4, Limit_reac_HB] / (-1 * v_1[Limit_reac_HB])
)

HB_Reactor_MB = Equation(
    container=m,
    name="HB_Reactor_MB",
    domain=[i],  # over components
    description="Mass balance over Haber-Bosch reactor"
)

HB_Reactor_MB[i] = F[5, i] == F[4, i] + v_1[i] * X_HB


# ===============================================================================#
#                                 || Flash  ||
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
NH3_Purity_Min = Parameter(
    container=m,
    name='Purity_Min',
    records=0.990,
    description="Minimum allowed liquid NH3 purity"
)

NH3_Purity_Max = Parameter(
    container=m,
    name='Purity_Max',
    records=0.996,
    description="Maximum allowed liquid NH3 purity"
)

# Purity Inequality Equations

NH3_Purity_Lower_Bound = Equation(
    container=m,
    name='Purity_Lower_Bound',
    description="Liquid NH3 purity must be >="
)
NH3_Purity_Lower_Bound[...] = F[6, 'NH3'] / Sum(i, F[6, i]) >= NH3_Purity_Min

NH3_Purity_Upper_Bound = Equation(
    container=m,
    name='Purity_Upper_Bound',
    description="Liquid NH3 purity must be <="
)
NH3_Purity_Upper_Bound[...] = F[6, 'NH3'] / Sum(i, F[6, i]) <= NH3_Purity_Max


# NH3 Recovery Definition
NH3_recovery = Variable(
    container=m,
    name="NH3_recovery",
    type="positive",
    description="NH3 recovery in flash unit"
)

NH3_Recov_Definition = Equation(
    container=m,
    name="NH3_Recov_Definition",
    description="Mass balance over Flash Separator"
)
NH3_Recov_Definition[...] = F[6, 'NH3'] == NH3_recovery * F[5, 'NH3']


vol_gases_only = Set(
    container=m,
    name='vol_gases_only',
    domain=[i],
    records=['H', 'Ar'],
    description="Permanent gases forced to vapor phase"
)

Zero_Gas_Flow = Equation(
    container=m,
    name='Zero_Gas_Flow',
    domain=[vol_gases_only],
    description="Permanent gases forced to vapor phase (F[6,i] == 0)"
)
Zero_Gas_Flow[vol_gases_only] = F[6, vol_gases_only] == 0

soluble_gases_only = Set(
    container=m,
    name='soluble_gases_only',
    domain=[i],
    records=['O', 'N', 'CO2'],
    description="Permanent gases forced to vapor phase"
)

Soluble_Gas_Flow = Equation(
    container=m,
    name='Soluble_Gas_Flow',
    domain=[soluble_gases_only],
    description="Gases soluble present in trace liquid phase"
)
Soluble_Gas_Flow[soluble_gases_only] = (
    F[6, soluble_gases_only] == (1 - 0.98) * F[5, soluble_gases_only]
)


H2O_recovery = Parameter(
    container=m,
    name="H2O_recovery",
    records=0.99,
    description="H2O recovery in flash unit"
)

H2O_Recov_Definition = Equation(
    container=m,
    name="H2O_Recov_Definition",
    description="H2O definition in flash unit"
)
H2O_Recov_Definition[...] = F[6, 'H2O'] == H2O_recovery * F[5, 'H2O']


# ===============================================================================#
#                                 || Split 1  ||
# ===============================================================================#

Split1MB = Equation(
    container=m,
    name="Split1MB",
    domain=[i],  # over components
    description="Mass balance over split 1"
)
Split1MB[i] = F[7, i] == F[8, i] + F[9, i]

HB_PurgeFraction = Parameter(
    container=m,
    name='HB_PurgeFraction',
    records=0.05,
    description="Purge fraction from HB recycle"
)


HB_PurgeFrac_def = Equation(
    container=m,
    domain=[i],
    name="HB_PurgeFrac_def",
    description="Definition of purge fraction from HB recycle"
)
HB_PurgeFrac_def[i] = F[8, i] == HB_PurgeFraction * F[7, i]


# ===============================================================================#
#                            || Maritime Transport ||
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
#                             || PRODUCTION TARGET ||
# ===============================================================================#

NH3_required_py = 46000  # ton/yr
NH3_MR = 17.031  # g/mol
hours_per_year = 8000  # hr/yr

NH3_mass_production = NH3_required_py * 1e3 / hours_per_year  # kg/yr
NH3_molar_production = NH3_mass_production / NH3_MR  # kmol/yr


# Parameter to hold the fixed production rate
NH3_Prod_Target = Parameter(
    container=m,
    name='NH3_Prod_Target',
    records=NH3_molar_production,
    description="Target NH3 production rate [kmol/hr] in Stream 6"
)

# Equation to fix the molar flow of NH3 in stream 6
Prod_Target_Constraint = Equation(
    container=m,
    name='Prod_Target_Constraint',
    description="Fixes NH3 output flowrate in Stream 6"
)

# F[6, 'NH3'] must equal the target amount.
Prod_Target_Constraint[...] = F[6, 'NH3'] == NH3_Prod_Target


# ===============================================================================#
#                             || INITIALIZATION FIX ||
# ===============================================================================#

# 1. Initialize all flows F to a small non-zero value for safety
F.l[j, i] = 1e-6  # Small non-zero initial guess for all flows


# 2. Fix flows of components that are known to be zero
# (e.g., NH3 in streams before production)

fix_values(F[2, 'NH3'], 0.0)
fix_values(F[1, 'NH3'], 0.0)

# Set initial recovery variables
NH3_recovery.lo[...] = 0.8
NH3_recovery.up[...] = 0.9

# H2O_recovery.lo[...] = 1e-6
# H2O_recovery.up[...] = 1.0


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

ObjFunc[...] = z == Sum(i, F[1, i] + F[2, i])  # Minimize fresh feed


# Define the Model
HB_Process_Model = Model(
    container=m,
    name='Haber_Bosch_Process',

    # 1. Set the objective to the H2 Fresh Feed flow
    objective=z,
    sense=Sense.MIN,

    # 2. List all required equations
    equations=m.getEquations(),

    # Use NLP because of the non-linear purity
    problem=Problem.NLP
)

# Solve the Model
print(HB_Process_Model.solve())
# Display results

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
