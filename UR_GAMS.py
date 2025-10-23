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
        'CO2',
        'AC',
        'UREA',
        'H2O'
    ],
    description="Involved chemical components"
)

# No. of streams in process
N_streams_start = 10
N_streams_end = 22

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
#                        || Urea Synthesis Loop (10+11+16->12+13)||
# ===============================================================================#
reactions = Set(
    container=m,
    name='reactions',
    records=[1, 2],
    description="1: Carbamate Formation, 2: Urea Synthesis"
)

# Stoichiometric Coefficients (Example for Carbamate/Urea)
# R1: CO2 + 2NH3 <-> AC (Ammonium Carbamate)
# R2: AC <-> UREA + H2O
nu = Parameter(
    container=m,
    domain=[i, reactions],
    name='nu',
    records=[
        # --- Reaction 1 (Carbamate Formation) ---
        ('CO2', 1, -1.0),  # CO2 consumed
        ('NH3', 1, -2.0),  # NH3 consumed
        ('AC', 1, 1.0),   # AC produced
        # The remaining components (UREA, H2O) are 0 in R1 by omission.

        # --- Reaction 2 (Urea Synthesis) ---
        ('AC', 2, -1.0),  # AC consumed
        ('UREA', 2, 1.0),  # UREA produced
        ('H2O', 2, 1.0)   # H2O produced
        # The remaining components (CO2, NH3) are 0 in R2 by omission.
    ],
    description="Stoichiometric coefficients for R1 and R2"
)

# Variables for reaction extent
extent = Variable(
    container=m,
    name="extent_of_reaction",
    domain=reactions,
    type="free",  # Can be negative for equilibrium/reversible reactions
    description="Extent of reaction for each reaction"
)

# 1. Component Mole Balances (i MBs + Extent)
USLmb = Equation(
    domain=i,
    container=m,
    name="USLmb",
    description="Urea Synthesis Loop Material Balance"
)
USLmb[i] = F['10', i] + F['11', i] + F['16', i] + \
           Sum(reactions, nu[i, reactions] * extent[reactions]) \
           == F['12', i] + F['13', i]


# 2. CO2 Feed Specification (5 Co2 Feed Spec)
S11_NoCO2 = Set(
    container=m,
    domain=i,
    name='S11_NoCO2',
    records=['AC', 'UREA', 'H2O', 'NH3'],
    description="Components not in the CO2 feed"
)
S11_CO2_Spec = Equation(
    container=m,
    domain=S11_NoCO2,  # i excluding CO2 and NH3
    name="S11_CO2_Spec",
    description="Stream 11 is pure CO2 feed"
)
# F_11,i = 0 if i != CO2, NH3
S11_CO2_Spec[S11_NoCO2] = F['11', S11_NoCO2] == 0.0  


# 3. SETgases Only in stream 12 (4 SETgases Only in stream 12)
SET_NoCO2NH3 = Set(
    container=m,
    domain=i,
    name='SET_NoCO2NH3',
    records=['AC', 'UREA', 'H2O'],
    description="Components not CO2, NH3"
)
S12_SETgases_Spec = Equation(
    container=m,
    domain=SET_NoCO2NH3,  # i excluding CO2 and NH3
    name="S12_SETgases_Spec",
    description="Only SETgases (CO2, NH3) leave in stream 12 purge"
)
# F_12,i = 0 if i != CO2, NH3
S12_SETgases_Spec[SET_NoCO2NH3] = F['12', SET_NoCO2NH3] == 0.0


# 4. Recycle Comp (4 Recycle Comp)
SET_NoACNH3CO2 = Set(
    container=m,
    domain=i,
    name='SET_NoACNH3CO2',
    records=['UREA', 'H2O'],
    description="Components not recycled (not AC, NH3)"
)
S16_Comp_Spec = Equation(
    container=m,
    domain=SET_NoACNH3CO2,  # i excluding AC and NH3
    name="S16_Comp_Spec",
    description="Only AC and NH3 are recycled in stream 16"
)
# F_16,i = 0 if i != AC, NH3
S16_Comp_Spec[SET_NoACNH3CO2] = F['16', SET_NoACNH3CO2] == 0.0


# 5. Purge Recovery (2 recoveries)
SETgases = Set(
    container=m,
    domain=i,
    name='SETgases',
    records=['CO2', 'NH3'],
    description="Gaseous components leaving in stream 12"
)
X_SynPurge = Parameter(
    container=m,
    name='X_SynPurge',
    records=0.05,  # Example value
    description="Fraction of (CO2+NH3) purged from stream 12"
)
PurgeRecovery = Equation(
    container=m,
    domain=SETgases,  # NH3 and CO2
    name="PurgeRecovery",
    description="Fraction of CO2 and NH3 purged in stream 12"
)
# F_12, (CO2/NH3) / (F_12,(CO2/NH3) + F_13,(CO2/NH3)) = X_SynPurge
# Set for components expected in SETgases (stream 12)
PurgeRecovery[SETgases] = (
    F['12', SETgases] == X_SynPurge * (F['12', SETgases] + F['13', SETgases])
)

# === REVISED EXTENT DEFINITION (Correct for Series Reactor) ===
OC_CO2 = Parameter(m, name='OC_CO2', records=0.75,
                   description="Target fractional conversion of CO2 (for R1)")
OC_AC = Parameter(m, name='OC_AC', records=0.60,
                  description="Target fractional conversion of AC (for R2)")
# R1: Extent based on total CO2 consumed
ExtentDef_R1 = Equation(
    container=m,
    name="ExtentDef_R1",
    description=(
        "xi_1 based on fractional conversion of total CO2 fed (F11+F16)"
    )
)
# xi_1 = OC_CO2 * (F_11,CO2 + F_16,CO2)
ExtentDef_R1[...] = extent['1'] == OC_CO2 * (F['11', 'CO2'] + F['16', 'CO2'])


# R2: Extent based on fractional conversion of total AC available
ExtentDef_R2 = Equation(
    container=m,
    name="ExtentDef_R2",
    description=(
        "xi_2 based on fractional conversion of AC available "
        "(F16_AC + xi_1)"
    )
)
# xi_2 = OC_AC * (F_16,AC + xi_1)
ExtentDef_R2[...] = extent['2'] == OC_AC * (F['16', 'AC'] + extent['1'])

# NEW - NH3:CO2 Molar Ratio in Feed Specification
R_NH3_CO2_Feed = 3 # Molar ratio of NH3 to CO2 in feeds to loop

NH3_CO2_FeedSpec = Equation(
    container=m,
    name="NH3_CO2_FeedSpec",
    description="NH3 to CO2 molar ratio in feeds to Urea Synthesis Loop"
)
NH3_CO2_FeedSpec[...] = R_NH3_CO2_Feed*(F['11', 'CO2']+F['16', 'CO2']) \
                        == F['10', 'NH3'] + F['16', 'NH3']

# ===============================================================================#
#                             || Decomposers (13->14+15) ||
# ===============================================================================#

# --- Reaction and Extent (Scalar Definitions) ---
# R: AC -> 2NH3 + CO2
nu_decomp = Parameter(
    container=m,
    domain=i,
    name='nu_decomp',
    records=[
        ('AC', -1.0),
        ('NH3', 2.0),
        ('CO2', 1.0)
    ],
    description="Stoichiometric coefficients for AC decomposition"
)

extent_decomp = Variable(
    container=m,
    name="extent_decomp",
    type="positive",
    description="Scalar extent of AC decomposition"
)

# --- Recovery Parameters ---
r_decomp_SETgases = Parameter(
    container=m,
    domain=SETgases,
    name='r_decomp_SETgases',
    records=[('CO2', 0.99), ('NH3', 0.99)],  # Example values
    description="Fraction of CO2 and NH3 recovered in stream 15"
)

# --- Equations ---

# 1. Component Mole Balances (i MBs + Scalar Extent)
DecomposersMB = Equation(
    domain=i,
    container=m,
    name="DecomposersMB",
    description="Decomposers Unit Material Balance"
)
DecomposersMB[i] = (
    F['13', i]
    + nu_decomp[i] * extent_decomp  # Scalar extent
    == F['14', i] + F['15', i]
)


# 2. Recoveries for CO2 and NH3 (2 equations)
DecompRecovery = Equation(
    container=m,
    domain=SETgases,
    name="DecompRecovery",
    description="Fraction of CO2 and NH3 recovered in stream 15"
)
DecompRecovery[SETgases] = (
    F['15', SETgases] == r_decomp_SETgases[SETgases]
    * (F['15', SETgases] + F['14', SETgases])
)


# 3. Urea & Water Split (UREA and H2O do not leave in overhead stream 15)
SET_UREAH2O = Set(
    container=m,
    domain=i,
    name='SET_UREAH2O',
    records=['UREA', 'H2O'],
    description="Components (UREA, H2O) that should not leave in stream 15"
)
S15_NoUreaH2O = Equation(
    container=m,
    domain=SET_UREAH2O,
    name="S15_NoUreaH2O",
    description="UREA and H2O go to stream 14 (bottoms)"
)
S15_NoUreaH2O[SET_UREAH2O] = F['15', SET_UREAH2O] == 0.0


# 4. Conversion (All feed AC Converted, limiting reactant)
DecompConversion = Equation(
    container=m,
    name="DecompConversion",
    description="All feed AC converted (AC is limiting reactant)"
)
DecompConversion[...] = extent_decomp == 1*F['13', 'AC']


# 5. Fix AC to ensure numerical stability
FixACSplit = Equation(
    container=m,
    name="FixACSplit",
    description="Fix F_14,AC = 0 to prevent numerical instability in splits"
)
FixACSplit[...] = F['14', 'AC'] == 0.0

# ===============================================================================#
#                           || Concentrators (14->17+18) ||
# ===============================================================================#

ConcentratorsMB = Equation(
    domain=i,
    container=m,
    name="ConcentratorsMB",
    description="Concentrators Unit Material Balance"
)
ConcentratorsMB[i] = (
    F['14', i]
    == F['17', i] + F['18', i]
)

r_concentrators = Parameter(
    container=m,
    domain=SET_UREAH2O,
    name='r_concentrators',
    records=[('H2O', 0.01), ('UREA', 0.99)],  # Example values
    description="Fraction of H2O and Urea recovered in stream 17"
)
ConcentratorRecovery = Equation(
    container=m,
    domain=SET_UREAH2O,
    name="ConcentratorRecovery",
    description="Fraction of H2O and Urea recovered in stream 17"
)
ConcentratorRecovery[SET_UREAH2O] = (
    F['17', SET_UREAH2O] == r_concentrators[SET_UREAH2O]
    * (F['14', SET_UREAH2O])
)

SET_NH3CO2 = Set(
    container=m,
    domain=i,
    name='SET_NH3CO2',
    records=['NH3', 'CO2'],
    description="Components (NH3, CO2) that should not leave in stream 17"
)

S17_NoNH3CO2 = Equation(
    container=m,
    domain=SET_NH3CO2,
    name="S17_NoNH3CO2",
    description="NH3 and CO2 go to stream 18 (bottoms)"
)
S17_NoNH3CO2[SET_NH3CO2] = F['14', SET_NH3CO2] == F['18', SET_NH3CO2]

# ===============================================================================#
#                           || Granulator (17->21+22) ||
# ===============================================================================#

X_GranLoss = Parameter(
    container=m,
    name='X_GranLoss',
    records=0.01,  # Example value
    description="Fractional loss of components in granulation"
)

GranulatorMB = Equation(
    domain=i,
    container=m,
    name="GranulatorMB",
    description="Granulator Unit Material Balance"
)
GranulatorMB[i] = (
    F['17', i]
    == F['21', i] + F['22', i]
)
GranulatorLoss = Equation(
    container=m,
    domain=i,
    name="GranulatorLoss",
    description="Fractional loss of components in granulation"
)
GranulatorLoss[i] = (
    F['22', i] == X_GranLoss * F['17', i]
)

# ===============================================================================#
#                              || PCT (18->19+20) ||
# ===============================================================================#

# R: UREA + H2O -> 2NH3 + CO2
nu_hydrolysis = Parameter(
    container=m,
    domain=i,
    name='nu_hydrolysis',
    records=[
        ('UREA', -1.0),
        ('H2O', -1.0),
        ('NH3', 2.0),
        ('CO2', 1.0)
    ]
)

# --- Reaction and Extent (Scalar Definitions) ---
extent_hydrolysis = Variable(
    container=m,
    name="extent_hydrolysis",
    type="positive",
    description="Scalar extent of Urea Hydrolysis in PCT"
)
extent_hydrolysis_def = Equation(
    container=m,
    name="extent_hydrolysis_def",
    description="Definition of extent of urea hydrolysis"
)
extent_hydrolysis_def[...] = extent_hydrolysis == F['18', 'UREA']

# --- Component Balances ---
PCT_MB = Equation(
    domain=i,
    container=m,
    name="PCT_MB",
    description="PCT Unit Material Balance"
)
PCT_MB[i] = (
    F['18', i] + nu_hydrolysis[i] * extent_hydrolysis
    == F['19', i] + F['20', i]
)

# 2 recovery params for NH3 and CO2 in stream 20
r_PCT_SETgases = Parameter(
    container=m,
    domain=SETgases,
    name='r_PCT_SETgases',
    records=[('CO2', 0.99), ('NH3', 0.99)],  # Example values
    description="Fraction of CO2 and NH3 recovered in stream 20"
)

S20_Gas_Recovery = Equation(
    container=m,
    domain=SETgases,
    name="S20_Gas_Recovery",
    description="Fraction of CO2 and NH3 recovered in stream 20"
)
S20_Gas_Recovery[SETgases] = (
    F['20', SETgases] == r_PCT_SETgases[SETgases]
    * (F['19', SETgases] + F['20', SETgases])
)

S19_NoUREAH2O = Equation(
    container=m,
    domain=SET_UREAH2O,
    name="S19_NoUREAH2O",
    description="UREA and H2O go to stream 19 (bottoms)"
)
S19_NoUREAH2O[SET_UREAH2O] = F['19', SET_UREAH2O] == 0.0
# ===============================================================================#
#                           || Absorbers (15+19->16) ||
# ===============================================================================#
extent_absorb = Variable(
    container=m,
    name="extent_absorb",
    type="positive",
    description="Scalar extent of Absorption"
)

absorber_efficiency = 0.90
extent_absorb_def = Equation(
    container=m,
    name="extent_absorb_def",
    description="Definition of extent of absorption"
)
extent_absorb_def[...] = extent_absorb == \
                         absorber_efficiency*(F['19', 'CO2'] + F['15', 'CO2'])

AbsorbersMB = Equation(
    domain=i,
    container=m,
    name="AbsorbersMB",
    description="Absorbers Unit Material Balance"
)
AbsorbersMB[i] = (
    F['15', i] + F['19', i] + nu[i, 1] * extent_absorb
    == F['16', i]
)


# ===============================================================================#
#                             || INLET REQUIREMENTS ||
# ===============================================================================#
# ===============================================================================#
#                             || INITIALIZATION FIX ||
# ===============================================================================#
# 1. Initialize all flows F to a small non-zero value for safety
F.l[j, i] = 1e-6  # Small non-zero initial guess for all flows

NH3_input = 135.015  # kmol/hr
fix_values(F[10, 'NH3'], NH3_input)
fix_values(F[10, 'H2O'], 0.0)
fix_values(F[10, 'UREA'], 0.0)
fix_values(F[10, 'AC'], 0.0)
fix_values(F[10, 'CO2'], 0.0)

extent.l['1'] = 1
extent.l['2'] = 1
extent_decomp.l[...] = 1
extent_hydrolysis.l[...] = 1
extent_absorb.l[...] = 1


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
ObjFunc[...] = z == Sum(i, F['21', i])  # Maximize total urea production

# Define the Model
UR_Prodution_Model = Model(
    container=m,
    name='UR_Prodution_Model',

    # 1. Set the objective to maximize urea production
    objective=z,
    sense=Sense.MAX,

    # 2. List all required equations
    equations=m.getEquations(),

    # Use NLP because of the non-linear purity
    problem=Problem.NLP
)

print(UR_Prodution_Model.solve())

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
