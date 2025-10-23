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
    records=list(range(1, N_streams_end + 1)),
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


# 3. Gases Only in stream 12 (4 Gases Only in stream 12)
S11_NoCO2NH3 = Set(
    container=m,
    domain=i,
    name='S11_NoCO2NH3',
    records=['AC', 'UREA', 'H2O'],
    description="Components not CO2, NH3"
)
S12_Gases_Spec = Equation(
    container=m,
    domain=S11_NoCO2NH3,  # i excluding CO2 and NH3
    name="S12_Gases_Spec",
    description="Only gases (CO2, NH3) leave in stream 12 purge"
)
# F_12,i = 0 if i != CO2, NH3
S12_Gases_Spec[S11_NoCO2NH3] = F['12', S11_NoCO2NH3] == 0.0


# 4. Recycle Comp (4 Recycle Comp)
S16_NoACNH3 = Set(
    container=m,
    domain=i,
    name='S16_NoACNH3',
    records=['CO2', 'UREA', 'H2O'],
    description="Components not recycled (not AC, NH3)"
)
S16_Comp_Spec = Equation(
    container=m,
    domain=S16_NoACNH3,  # i excluding AC and NH3
    name="S16_Comp_Spec",
    description="Only AC and NH3 are recycled in stream 16"
)
# F_16,i = 0 if i != AC, NH3
S16_Comp_Spec[S16_NoACNH3] = F['16', S16_NoACNH3] == 0.0


# 5. Purge Recovery (2 recoveries)
gases = Set(
    container=m,
    domain=i,
    name='gases',
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
    domain=gases,  # NH3 and CO2
    name="PurgeRecovery",
    description="Fraction of CO2 and NH3 purged in stream 12"
)
# F_12, (CO2/NH3) / (F_12,(CO2/NH3) + F_13,(CO2/NH3)) = X_SynPurge
# Set for components expected in gases (stream 12)
PurgeRecovery[gases] = (
    F['12', gases] == X_SynPurge * (F['12', gases] + F['13', gases])
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
r_decomp_gases = Parameter(
    container=m,
    domain=gases,
    name='r_decomp_gases',
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
    domain=gases,
    name="DecompRecovery",
    description="Fraction of CO2 and NH3 recovered in stream 15"
)
DecompRecovery[gases] = (
    F['15', gases] == r_decomp_gases[gases]
    * (F['15', gases] + F['14', gases])
)


# 3. Urea & Water Split (UREA and H2O do not leave in overhead stream 15)
S15_UREAH2O = Set(
    container=m,
    domain=i,
    name='S15_UREAH2O',
    records=['UREA', 'H2O'],
    description="Components (UREA, H2O) that should not leave in stream 15"
)
S15_NoUreaH2O = Equation(
    container=m,
    domain=S15_UREAH2O,
    name="S15_NoUreaH2O",
    description="UREA and H2O go to stream 14 (bottoms)"
)
S15_NoUreaH2O[i] = F['14', i] == 0.0


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
