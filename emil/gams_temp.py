from gamspy import (
    Container,
    Set,
    Parameter,
    Variable,
    Equation,
    Model,
    Sense,
    Card,
    Domain,
    Sum,
    SetElement,
    Constant,
)
import numpy as np
import pandas as pd
from io import StringIO

# Initialize the GAMSpy Container
c = Container()

## ------------------------------ SETS -------------------------------

# *Initialisation of sets for components of the process, c, streams, s and Antoine's coefficients, coeff.
# *Here, dme - dimethyl ether, m - methanol, co - carbon monoxide, h - hydrogen, mea - methyl acetate.
comp = Set(
    c,
    domain="*",
    records=[
        "dme",
        "m",
        "co",
        "h",
        "mea",
    ],
    description="Components",
)

# *Stream 1 - DME feed, 2 - CO feed, 3 - mixer outlet, 4 - reactor effluent, 5 - flash drum vapour output,
# *6 - flash drum liquid output, 7 - gas recycle, 8 - purge, 9 - distillate vapour, 10 - product stream.
stream = Set(
    s,
    domain="*",
    records=[str(i) for i in range(1, 11)],
    description="Streams",
)

coeff = Set(
    coeff, domain="*", records=["A", "B", "C"], description="Antoine's coefficients"
)

## ------------------------------ VARIABLES & PARAMETERS -------------------------------

# *Initialisation of a variable representing individual flow rate of component c in the stream s.
f = Variable(
    f, domain=[stream, comp], type="Positive", description="Flow rate of component c in stream s"
)

# *Initialisation of the value for total CO feed stream flow rate.
totalCO = Parameter(totalCO, description="Total CO feed stream flow rate")
totalCO.set_val(268)

# *Initialisation of pressure and conversion of a reactor.
# *Nu(c) - stoichiometric coefficient is given for each of the reactants
P = Parameter(P, description="Reactor Pressure")
P.set_val(30)
conv = Parameter(conv, description="Reactor Conversion")
# conv = 0.013373*P + 0.040682 is calculated inside the loop

# *Xi - extent of reaction is set as a variable.
xi = Variable(xi, type="Positive", description="Extent of reaction")

# *Definition of purity and purity equation to remove an unknown in the product stream.
purity = Parameter(purity, description="Product purity")
purity.set_val(0.99)

# *Initialisation of Flash Drum Temperature and Pressure, as well as Vapour Pressure for every component.
FlashP = Parameter(FlashP, description="Flash Drum Pressure")
FlashP.set_val(29)

# pstar1 is a parameter but its value depends on FlashT, so we define FlashT first.
FlashT = Variable(FlashT, description="Flash Drum Temperature")
FlashT.set_bounds(308.15, 320.15)
FlashT.set_starting_point(308.15)

pstar1 = Parameter(pstar1, domain=comp, description="Vapour Pressure at Flash Drum T")

# *Initialisation of variables for the flash drum calculations,
x = Variable(x, domain=comp, description="Liquid composition")
y = Variable(y, domain=comp, description="Gas composition")
z = Variable(z, domain=comp, description="Molar composition")
L = Variable(L, description="Total liquid flow")
V = Variable(V, description="Total Vapour flow")
Effl = Variable(Effl, description="Effluent flow rate (Stream 4)")

# *Initialisation of parameters such as distillation temperature and relative volatility
DTemp = Parameter(DTemp, description="Distillation Temperature")
DTemp.set_val(308.15)
pstar2 = Parameter(pstar2, domain=comp, description="Vapour Pressure at Distillation T")
Al = Parameter(Al, description="Relative Volatility (DME/MEA)")

# *Initialisation of total distillate flow rate and recovery.
totdist = Variable(totdist, description="Total distillate flow rate (Stream 9)")
r = Variable(r, description="Recovery")
r.set_bounds(0.99, 0.998)

# *Initialisation of the spit fraction
sf = Parameter(sf, description="Split fraction (Stream 8 / Stream 5)")
sf.set_val(0.05)

# *Definition of operating annual hours, cost of the feed specified in the brief and prices for CO and MEA.
hpy = Parameter(hpy, description="Hours per year")
hpy.set_val(8150)
costfeed = Parameter(costfeed, description="Cost of Feed 1 (DME/Methanol mixture)")
costfeed.set_val(78.2)
price = Parameter(price, domain=comp, description="Prices/Costs for components")
price.setRecords(pd.DataFrame([["co", 13.2], ["mea", 135.0]], columns=["c", "value"]))

# *Initialisation of net value of raw materials cost compared to product (MEA) selling price.
netv = Variable(netv, description="Net value of raw materials cost vs product price (MM$/yr)")

# *Initialisation of capital costs, compressor operation costs and utility costs.
ReactorCost = Variable(ReactorCost, description="Reactor cost (MM$)")
FlashCost = Variable(FlashCost, description="Flash Vessel cost (MM$)")
DistilCost = Variable(DistilCost, description="Distillation unit cost (MM$)")

Comp1Cost = Variable(Comp1Cost, description="Compressor 1 cost (MM$) - CO Feed")
Comp2Cost = Variable(Comp2Cost, description="Compressor 2 cost (MM$) - Recycle Gas")
Comp3Cost = Variable(Comp3Cost, description="Compressor 3 cost (MM$) - DME Feed")
Comp4Cost = Variable(Comp4Cost, description="Compressor 4 cost (MM$) - Distillate")

Comp1WCost = Variable(Comp1WCost, description="Compressor 1 work cost (MW)")
Comp2WCost = Variable(Comp2WCost, description="Compressor 2 work cost (MW)")
Comp3WCost = Variable(Comp3WCost, description="Compressor 3 work cost (MW)")
Comp4WCost = Variable(Comp4WCost, description="Compressor 4 work cost (MW)")

UtilityCost = Variable(UtilityCost, description="Utility Costs (MM$)")

TotFlRec = Parameter(TotFlRec, description="Total flow of recycle stream (Initial guess)")

# *Initialisation of the unit costs of units which is simply a sum of the individual TACs.
unitc = Variable(unitc, description="Total Annualized Capital Cost (MM$)")

# *Initialisation of the economic potential for the process and its equation.
econpot = Variable(econpot, description="Economic Potential (MM$/yr)")
econpot.set_upper_bound(1.0e10)

# *Reformulation of Economic potential into an objective function which is to be maximised.
obj = Variable(obj, description="Objective function = Economic Potential")

## ------------------------------ DATA INITIALISATION -------------------------------

# nu(c) - stoichiometric coefficient
nu = Parameter(
    nu,
    domain=comp,
    records=pd.DataFrame(
        [
            ["dme", -1],
            ["co", -1],
            ["mea", 1],
            ["m", 0],
            ["h", 0],
        ],
        columns=["c", "value"],
    ),
    description="Stoichiometric coefficient",
)

# Antoine's coefficients
antoine_data = """
      A        B          C
h     3.54314  99.395     7.726
co    3.36515  230.272   -13.15
dme   4.11475  894.669   -30.604
mea   4.20364  1164.426  -52.69
m     5.15853  1569.613  -34.846
"""
antoine_df = pd.read_csv(StringIO(antoine_data), sep=r"\s+", index_col=0)
antoine = Parameter(antoine, domain=[comp, coeff], description="Antoine's coefficients")
for c_name in antoine_df.index:
    for coeff_name in antoine_df.columns:
        if c_name in comp.get_domain_elements() and coeff_name in coeff.get_domain_elements():
            antoine.set_value([c_name, coeff_name], antoine_df.loc[c_name, coeff_name])

## ------------------------------ FIXED VALUES & GUESSES -------------------------------

# Setting fixed values for feed 1
f.set_value(["1", "dme"], 0.999 * 260)
f.set_value(["1", "m"], 0.001 * 260)
f.set_value(["1", "co"], 0)
f.set_value(["1", "h"], 0)
f.set_value(["1", "mea"], 0)

# Setting lower bounds and fixed values for feed 2
f.set_lower_bound(["2", "co"], 0.98 * 260)
f.set_lower_bound(["2", "h"], 0.02 * 260)
f.set_value(["2", "dme"], 0)
f.set_value(["2", "m"], 0)
f.set_value(["2", "mea"], 0)

# Guesses for CO feed stream component flow rates
f.set_starting_point(["2", "co"], 281)
f.set_starting_point(["2", "h"], 5.62)

# Flash Drum Vapour Pressure calculation (using initial guess for FlashT.l)
# pstar1(c) = 10**(antoine(c,'A') - antoine(c,'B') / (FlashT.l + antoine(c,'C')));
# We'll use the .l (level) value of FlashT for the parameter calculation before the solve
# The parameter will be updated inside the loop based on the previous solution's FlashT.l
def update_pstar1(flash_t_level):
    pstar1_data = {
        c_name: 10
        ** (
            antoine.get_value([c_name, "A"])
            - antoine.get_value([c_name, "B"])
            / (flash_t_level + antoine.get_value([c_name, "C"]))
        )
        for c_name in comp.get_domain_elements()
    }
    for c_name, val in pstar1_data.items():
        pstar1.set_value(c_name, val)

update_pstar1(FlashT.get_starting_point())

# Initial guesses for flash drum flow rates
# f.l('6',c) = f.l('4',c) * 0.5;
# f.l('5',c) = f.l('4',c) * 0.5;
# f.l('4',c) = f.l('3',c) + nu(c) * xi.l
# Since f('4',c) is an intermediate, we'll set f('5',c) and f('6',c) to an arbitrary small positive value for initial solution
f.set_starting_point(
    [stream.lead_by_name("5"), comp.lead()], 1.0
)  # Arbitrary small positive
f.set_starting_point(
    [stream.lead_by_name("6"), comp.lead()], 1.0
)  # Arbitrary small positive


# Distillation column vapour pressure and relative volatility
# pstar2(c) = 10**(antoine(c,'A') - antoine(c,'B') / (DTemp + antoine(c,'C')));
pstar2_data = {
    c_name: 10
    ** (
        antoine.get_value([c_name, "A"])
        - antoine.get_value([c_name, "B"])
        / (DTemp.get_value() + antoine.get_value([c_name, "C"]))
    )
    for c_name in comp.get_domain_elements()
}
for c_name, val in pstar2_data.items():
    pstar2.set_value(c_name, val)

# Al = (pstar2('dme'))/pstar2('mea');
Al.set_val(pstar2.get_value("dme") / pstar2.get_value("mea"))

# Fixed values for Distillation Column output
f.set_value(["10", "h"], 0)
f.set_value(["10", "co"], 0)
f.set_value(["9", "m"], 0)

# TotFlRec initial guess
TotFlRec.set_val(sum(f.get_starting_point([stream.lead_by_name("7"), comp.lead()])))

## ------------------------------ EQUATIONS -------------------------------

# Total CO feed flow rate definition
totalC = Equation(totalC, description="Total CO feed flow rate")
COfeedfr = Equation(COfeedfr, description="Fraction of CO in the feed")
RatioCO = Equation(RatioCO, description="Composition of H2 in the stream")

totalC.expr = totalCO == Sum(comp, f["2", comp])
COfeedfr.expr = f["2", "co"] == 0.98 * totalCO
RatioCO.expr = totalCO == 50 * f["2", "h"]

# Constraint of DME feed being greater than feed of CO
rfc = Equation(rfc, description="CO flow rate in stream 3 >= DME flow rate in stream 3")
rfc.expr = f["3", "co"] >= f["3", "dme"]

# Purity equation
purityP = Equation(purityP, description="Purity of MEA in product stream 10")
purityP.expr = f["10", "mea"] == purity * Sum(comp, f["10", comp])

# ------------------------------ MIXER -------------------------------

# Mixer mass balance equation.
mixermb = Equation(mixermb, domain=comp, description="Mixer mass balance")
mixermb.expr = (
    f["1", comp] + f["2", comp] + f["7", comp] + f["9", comp] == f["3", comp]
)

# ----------------------------- REACTOR ------------------------------

# Initialisation of equation for the conversion definition and reactor mass balance.
conversiondef = Equation(conversiondef, description="Conversion definition")
reactormb = Equation(reactormb, domain=comp, description="Reactor mass balance")

conversiondef.expr = xi == conv * f["3", "dme"]
reactormb.expr = f["4", comp] == f["3", comp] + nu[comp] * xi

# -------------------------------- FLASH DRUM --------------------------------

# Initialisation of equations using the variables set above, in addition to Raoult's law
# equation and flash drum mass balance.
FlashInput = Equation(FlashInput, description="Flash Input Flow")
LEQ = Equation(LEQ, description="Total Liquid Flow")
VEQ = Equation(VEQ, description="Total Vapour Flow")

xeq = Equation(xeq, domain=comp, description="Liquid flow definition")
yeq = Equation(yeq, domain=comp, description="Vapour flow definition")
zeq = Equation(zeq, domain=comp, description="Effluent composition definition")

zsum = Equation(zsum, description="Sum of molar composition z is 1")
xsum = Equation(xsum, description="Sum of liquid composition x is 1")
ysum = Equation(ysum, description="Sum of gas composition y is 1")

Raoults = Equation(Raoults, domain=comp, description="Raoult's Law")
Flashb = Equation(Flashb, domain=comp, description="Flash drum mass balance")

FlashInput.expr = Effl == Sum(comp, f["4", comp])
LEQ.expr = L == Sum(comp, f["6", comp])
VEQ.expr = V == Sum(comp, f["5", comp])

xeq.expr = f["6", comp] == x[comp] * L
yeq.expr = f["5", comp] == y[comp] * V
zeq.expr = f["4", comp] == z[comp] * Effl

zsum.expr = Sum(comp, z[comp]) == 1
xsum.expr = Sum(comp, x[comp]) == 1
ysum.expr = Sum(comp, y[comp]) == 1

Raoults.expr = y[comp] * FlashP == x[comp] * pstar1[comp]
Flashb.expr = z[comp] * Effl == x[comp] * L + y[comp] * V

# -------------------------------- DISTILLATION COLUMN --------------------------------

# Initialisation of the distillation column mass balance, recovery equations and total distillate.
distilmb = Equation(distilmb, domain=comp, description="Distillation mass balance")
LK = Equation(LK, description="Light Key (DME) Recovery")
HK = Equation(HK, description="Heavy Key (MEA) Recovery")
distillate = Equation(distillate, description="Total Distillate Flow")

distilmb.expr = f["6", comp] == f["9", comp] + f["10", comp]
LK.expr = f["9", "dme"] == r * f["6", "dme"]
HK.expr = f["10", "mea"] == r * f["6", "mea"]
distillate.expr = totdist == Sum(comp, f["9", comp])

# ------------------------------------- SPLITTER --------------------------------------

# Initialisation of equations - mass balance and split fraction.
splitmb = Equation(splitmb, domain=comp, description="Splitter mass balance")
splitsf = Equation(splitsf, domain=comp, description="Split fraction definition")

splitmb.expr = f["5", comp] == f["7", comp] + f["8", comp]
splitsf.expr = f["8", comp] == sf * f["5", comp]

# ---------------------------------- ECONOMIC EVALUATION ----------------------------------

# Net value of raw materials cost compared to product (MEA) selling price.
materialscosts = Equation(
    materialscosts, description="Net value of raw materials vs product price"
)
materialscosts.expr = netv == (
    hpy
    * (
        (price["mea"] * Sum(comp, f["10", comp]))
        - (price["co"] * Sum(comp, f["2", comp]))
        - (costfeed * 260)
    )
) / 1000000

# Capital and Operating Costs
CostFlash = Equation(CostFlash, description="Flash Vessel cost")
CostDistil = Equation(CostDistil, description="Distillation unit cost")
CostReactor = Equation(CostReactor, description="Reactor cost")

CostComp1 = Equation(CostComp1, description="Compressor 1 cost (CO Feed)")
CostComp2 = Equation(CostComp2, description="Compressor 2 cost (Recycle Gas)")
CostComp3 = Equation(CostComp3, description="Compressor 3 cost (DME Feed)")
CostComp4 = Equation(CostComp4, description="Compressor 4 cost (Distillate)")

CostComp1W = Equation(CostComp1W, description="Compressor 1 work cost")
CostComp2W = Equation(CostComp2W, description="Compressor 2 work cost")
CostComp3W = Equation(CostComp3W, description="Compressor 3 work cost")
CostComp4W = Equation(CostComp4W, description="Compressor 4 work cost")

CostUtility = Equation(CostUtility, description="Utility Costs")

CostFlash.expr = FlashCost == 0.001 * Sum(comp, f["4", comp])
CostDistil.expr = (
    DistilCost * (100 * (1 - r)) * (Al - 1) == Sum(comp, f["6", comp]).sqrt()
)
CostReactor.expr = ReactorCost * 10**6 == 13248 * Sum(comp, f["3", comp])

# Note: The GAMS code for compressor costs (Comp1-4Cost) uses the sum of flow rates
# f('2',c), TotFlRec, f('1',c), f('9',c) respectively, which is implemented below:

CostComp1.expr = CostComp1 == 0.02 * Sum(comp, f["2", comp])
CostComp2.expr = CostComp2 == 0.02 * TotFlRec
CostComp3.expr = CostComp3 == 0.02 * Sum(comp, f["1", comp])
CostComp4.expr = CostComp4 == 0.02 * Sum(comp, f["9", comp])

# Compressor Work Costs
# Constants for the work equations:
# R = 8.314 J/mol.K
# gamma = 1.4, gamma/(gamma-1) = 1.4/0.4 = 3.5

# P_out / P_in Ratios (approximated based on GAMS code's final pressure and context):
# Comp1: P_in ~ 1 atm (or low), P_out ~ 33 bar. Ratio = 33 / P_in (GAMS uses 33**0.4/1.4 - 1, implying P_in ~ 1.01325)
# Comp2: P_in ~ 29 bar (FlashP), P_out ~ 33 bar. Ratio = 33/29
# Comp3: P_in ~ 1.01325 bar (Atmospheric), P_out ~ 33 bar. Ratio = 33/1.01325
# Comp4: P_in ~ 5.06625 bar, P_out ~ 33 bar. Ratio = 33/5.06625

# R_const = 8.314
# C_ratio = 1.4 / 0.4
# R_exp = 0.4 / 1.4

# Comp1: Flow ~ TotalCO, T=298.15K
CostComp1W.expr = CostComp1W == 1.562 * (
    (totalCO * (1 / 3600) * (1 / 1000))
    * (1.4 / 0.4)
    * (8.314)
    * (298.15)
    * ((33) ** (0.4 / 1.4) - 1)
)

# Comp2: Flow ~ 50 * Comp2Cost (Approximating Flow from Cost based on GAMS structure), T=308.15K
CostComp2W.expr = CostComp2W == 1.562 * (
    (50 * Comp2Cost * (1 / 3600) * (1 / 1000))
    * (1.4 / 0.4)
    * (8.314)
    * (308.15)
    * ((33 / 29) ** (0.4 / 1.4) - 1)
)

# Comp3: Flow ~ 50 * Comp3Cost (Approximating Flow from Cost), T=308.15K
CostComp3W.expr = CostComp3W == 1.562 * (
    (50 * Comp3Cost * (1 / 3600) * (1 / 1000))
    * (1.4 / 0.4)
    * (8.314)
    * (308.15)
    * ((33 / 1.01325) ** (0.4 / 1.4) - 1)
)

# Comp4: Flow ~ 50 * Comp4Cost (Approximating Flow from Cost), T=308.15K
CostComp4W.expr = CostComp4W == 1.562 * (
    (50 * Comp4Cost * (1 / 3600) * (1 / 1000))
    * (1.4 / 0.4)
    * (8.314)
    * (308.15)
    * ((50 * Comp4Cost * (1 / 3600) * (1 / 1000))
    * (1.4 / 0.4)
    * (8.314)
    * (308.15)
    * ((33 / 5.06625) ** (0.4 / 1.4) - 1)
    )
)

# Utility Cost (Annual Cost of Compressor Work)
CostUtility.expr = (
    UtilityCost * 10**9
    == (Comp1WCost + Comp2WCost + Comp3WCost + Comp4WCost) * 3600 * 8150 * 17
)

# Total Annualized Capital Cost
costunits = Equation(costunits, description="Total Annualized Capital Cost")
costunits.expr = unitc == (
    ReactorCost
    + FlashCost
    + DistilCost
    + Comp1Cost
    + Comp2Cost
    + Comp3Cost
    + Comp4Cost
)

# Economic Potential
profit = Equation(profit, description="Economic Potential")
profit.expr = econpot == netv - unitc - UtilityCost

# Objective function
objective = Equation(objective, description="Objective function = Economic Potential")
objective.expr = obj == econpot

# Model definition
process = Model(
    c,
    name="process",
    equations=c.get<ctrl61>_all_equations(),
    sense=Sense.MAX,
    objective=obj,
    problem="nlp",
)

## ------------------------------ LOOPING & SOLVE -------------------------------

# Prepare data structure for results
results_columns = [
    "#",
    "totalCO",
    "totalDist",
    "sf",
    "Pr",
    "conv",
    "r",
    "obj",
    "FlashT",
    "status",
]
results_df = pd.DataFrame(columns=results_columns)
w = 0

# Convert the GAMS loops into nested Python loops
totalCO_values = np.arange(260, 290 + 2, 2)  # 260 to 290 by 2
sf_values = np.arange(0.0, 0.2 + 0.02, 0.02)  # 0.0 to 0.2 by 0.02
P_values = np.arange(20, 30 + 1, 1)  # 20 to 30 by 1

print(
    "{:^12} {:^10} {:^10} {:^10} {:^10} {:^10} {:^10} {:^12} {:^10} {:^10}".format(
        *results_columns
    )
)
print("-" * 116)

for tc_val in totalCO_values:
    totalCO.set_val(tc_val)

    for s_val in sf_values:
        sf.set_val(s_val)

        for p_val in P_values:
            P.set_val(p_val)

            # Recalculation of conversion
            conv_val = 0.013373 * p_val + 0.040682
            conv.set_val(conv_val)

            # Setting level values for all stream flow rates.
            # f.l(s ,c) = 260;
            f.set_starting_point([stream.lead(), comp.lead()], 260)
            # f.l('7',c) = sf * 260;
            f.set_starting_point([stream.lead_by_name("7"), comp.lead()], s_val * 260)

            # Update initial guess for recycle total flow
            TotFlRec.set_val(
                sum(f.get_starting_point([stream.lead_by_name("7"), comp.lead()]))
            )

            # Update the Flash Drum Vapour Pressure parameter pstar1 based on FlashT.l from the *previous* run's solution.
            # For the very first run, it uses the initial guess (308.15).
            update_pstar1(FlashT.get_starting_point())

            # Updating the counter.
            w += 1

            # Utilisation of NLP for non linear system and optimisation.
            process.solve()

            # Check model status
            # GAMS modelstat 2 is "Optimal"
            if (
                process.model_status in [2]
            ):
                econpot_l = econpot.get_latest_level()
                if econpot_l >= 0:
                    # Collect results
                    new_row = {
                        "#": w,
                        "totalCO": tc_val,
                        "totalDist": totdist.get_latest_level(),
                        "sf": s_val,
                        "Pr": p_val,
                        "conv": conv_val,
                        "r": r.get_latest_level(),
                        "obj": obj.get_latest_level(),
                        "FlashT": FlashT.get_latest_level(),
                        "status": process.model_status,
                    }

                    # Append to results DataFrame
                    results_df.loc[len(results_df)] = new_row
                    
                    # Print desired values (similar to GAMS 'put')
                    print(
                        "{:^12} {:^10.2f} {:^10.2f} {:^10.2f} {:^10.2f} {:^10.4f} {:^10.4f} {:^12.4f} {:^10.2f} {:^10}".format(
                            new_row["#"],
                            new_row["totalCO"],
                            new_row["totalDist"],
                            new_row["sf"],
                            new_row["Pr"],
                            new_row["conv"],
                            new_row["r"],
                            new_row["obj"],
                            new_row["FlashT"],
                            new_row["status"],
                        )
                    )

                    # Update starting points for the next iteration with the solved levels
                    f.set_starting_point(
                        f.get_latest_levels()
                    )
                    r.set_starting_point(r.get_latest_level())
                    FlashT.set_starting_point(FlashT.get_latest_level())
                    xi.set_starting_point(xi.get_latest_level())
                    
                    # TotFlRec needs to be updated for Comp2Cost for the next iteration
                    TotFlRec.set_val(
                        sum(f.get_latest_levels().loc[stream.lead_by_name("7"), :])
                    )


# Write results to file (similar to GAMS 'file results' and 'put results')
# The GAMS code writes to PDP-Results.data, we'll use that name.
results_df.to_csv("PDP-Results.data", sep="\t", index=False)

# Optional: Displaying the overall best result
if not results_df.empty:
    best_result = results_df.loc[results_df["obj"].idxmax()]
    print("\n" + "=" * 116)
    print("✨ Best Economic Potential Found (from successful, non-negative runs):")
    print(
        "{:^12} {:^10.2f} {:^10.2f} {:^10.2f} {:^10.2f} {:^10.4f} {:^10.4f} {:^12.4f} {:^10.2f} {:^10}".format(
            best_result["#"],
            best_result["totalCO"],
            best_result["totalDist"],
            best_result["sf"],
            best_result["Pr"],
            best_result["conv"],
            best_result["r"],
            best_result["obj"],
            best_result["FlashT"],
            best_result["status"],
        )
    )
    print("=" * 116)
else:
    print("\nNo successful, non-negative economic potential solutions were found.")