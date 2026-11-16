#!/usr/bin/env python3
"""Motor-case pressure check + combustion-pressure uncertainty + hydro-test ladder.

!! Geometry / material values below are EXAMPLE placeholders, NOT the team's design.
!! Replace the EXAMPLE block with the real drawing and mill-sheet / standard values.
Tags: [P] paper, [E] EXAMPLE placeholder (replace!), [A] assumption to verify, [D] derived.
Run: python3 calc/case_pressure.py
"""
import math, json, pathlib

OUT = pathlib.Path(__file__).parent

# ---------------- [E] EXAMPLE geometry / materials (REPLACE) ----------------
E = dict(
    Do=82.0,            # mm  outer dia  [P] Table 1 (paper's motor case)
    Di=70.0,            # mm  inner dia  [E] not given in the paper
    groove_depth=3.0,   # mm  O-ring groove depth cutting into the wall [E]
    sigma_y=240.0,      # MPa  6061-T6 yield, typical minimum  [A] verify against material sheet
    sigma_u=290.0,      # MPa  6061-T6 ultimate, typical minimum [A] verify
    D_seal=70.0,        # mm  diameter on which pressure acts on the closure [E]
    n_bolt=6, As_bolt=20.1, sigma_u_bolt=1200.0,  # M6, tensile stress area ~20.1 mm2, grade 12.9 [A][E]
    Pc_design=2.5,      # MPa [P] Table 2
    ratio_proof=1.5,    # whiteboard: "1.5倍?"
)

def lame_yield_pressure(ro, ri, sy):
    """Hoop stress at inner wall reaches sy (thick-wall Lame, conservative: uses hoop only)."""
    return sy * (ro**2 - ri**2) / (ro**2 + ri**2)

def barlow_burst(t, Do, su):
    """Thin-wall burst estimate P = 2 t su / Do."""
    return 2 * t * su / Do

ro, ri = E["Do"] / 2, E["Di"] / 2
t = ro - ri
modes = {}
modes["wall: first yield (Lame)"] = lame_yield_pressure(ro, ri, E["sigma_y"])
modes["wall: burst est. (Barlow, UTS)"] = barlow_burst(t, E["Do"], E["sigma_u"])
tg = t - E["groove_depth"]
modes["wall at O-ring groove: first yield (Lame)"] = lame_yield_pressure(ro, ri + E["groove_depth"], E["sigma_y"])
modes["wall at O-ring groove: burst est. (Barlow)"] = barlow_burst(tg, E["Do"], E["sigma_u"])
A_cl = math.pi / 4 * E["D_seal"] ** 2           # mm2
cap_bolt = E["n_bolt"] * E["As_bolt"] * E["sigma_u_bolt"]   # N  (tensile only, no thread strip, no preload)
modes["closure bolts: tensile break (no thread-strip check)"] = cap_bolt / A_cl
weakest = min(modes, key=modes.get)

Pproof = E["Pc_design"] * E["ratio_proof"]
res = dict(
    NOTE="EXAMPLE geometry. Not the team's design. Replace the [E] values.",
    inputs=E, wall_thickness_mm=t, wall_at_groove_mm=tg,
    closure_area_mm2=A_cl, closure_force_kN_at_Pc=E["Pc_design"] * A_cl / 1e3,
    failure_pressures_MPa=modes, weakest_mode=weakest,
    test_pressure_MPa=Pproof,
    margin_weakest_over_proof=modes[weakest] / Pproof,
    ladder_MPa=[round(x, 2) for x in (0.5 * Pproof, 0.75 * Pproof, Pproof)],
)

# ---------------- combustion-pressure uncertainty (whiteboard: "わからない") ----------------
# Pc = mdot * c* / At.  Baseline from calc/candy_motor.py: mdot=0.2133 kg/s, c*=1141 m/s (gamma-dependent!), At=97.4 mm2 -> 2.5 MPa
MDOT, CSTAR = 0.2133, 1141.0
rows = []
for d in (9.0, 10.0, 11.13, 12.0, 13.0):
    At = math.pi / 4 * (d * 1e-3) ** 2
    lo = 0.85 * MDOT * 0.85 * CSTAR / At / 1e6
    mid = MDOT * CSTAR / At / 1e6
    hi = 1.15 * MDOT * 1.15 * CSTAR / At / 1e6
    rows.append(dict(throat_mm=d, Pc_MPa_low=round(lo, 2), Pc_MPa_nominal=round(mid, 2), Pc_MPa_high=round(hi, 2)))
res["Pc_vs_throat(mdot,c* each +-15%; steady state, no ignition spike, no erosion)"] = rows
(OUT / "case_pressure_results.json").write_text(json.dumps(res, indent=2, ensure_ascii=False))
print(json.dumps(res, indent=2, ensure_ascii=False))
