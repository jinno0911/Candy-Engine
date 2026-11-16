#!/usr/bin/env python3
"""Candy hybrid rocket motor: design calculation based on
Wada et al., Trans. JSASS Aerospace Tech. Japan 16(6), 506-510 (2018).

Every number is tagged:  [P] = read from the paper (table/fig/page)
                         [A] = ASSUMPTION (not in the paper, needs verification)
                         [D] = derived here from the above
Run:  python3 calc/candy_motor.py   -> writes calc/results.json and calc/*.png
"""
import json, math, pathlib
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = pathlib.Path(__file__).parent
G0 = 9.80665            # standard gravity [m/s^2]  (definition)
PA = 101325.0           # sea-level pressure [Pa]  (standard)

# ---------------- [P] paper values ----------------
P = dict(
    Pc=2.5e6,            # Table 2: combustion pressure 2.5 MPa
    burn_time_plan=3.0,  # Table 2
    F_avg=351.8,         # Table 4 soft candy, average thrust [N]
    t_burn=2.62,         # Sec.2.3 text (Table 4 says 2.6 s)
    It_text=921.7,       # Sec.2.3 text: total impulse [N s]
    Isp=168.2,           # Table 4 soft candy [s]
    m_fuel=0.064,        # Sec.2.2: 64 g of fuel in chamber
    L_motor=0.339, D_out=0.082, D_in=0.050,   # Table 1
    G_ox_table2=70.0,    # Table 2 [kg/m2 s]
    G_ox_text=300.0,     # Sec.2.3 text [kg/m2 s]  <-- contradicts Table 2
    rocket_dry=10.0, rocket_D=0.150,          # Table 5
    alt_calc_candy=313.0, alt_meas_candy=319.0, t_apogee_fig7=9.0,  # Table5, Fig7 (read off ~9 s)
    A_wet=4.240, A_dry=3.670, A_D=0.089, A_alt_calc=1289.0,   # Table 6
    B_wet=4.652, B_dry=4.322, B_D=0.089, B_alt_calc=1187.0,   # Table 7
)
# ---------------- [A] assumptions ----------------
A = dict(
    gamma=1.2,           # [A] typical exhaust specific-heat ratio, NOT from paper
    rho_N2O_liq=786.0,   # [A] liquid N2O density ~20 C [kg/m3]; verify with a property table
    tank_margin=1.20,    # [A] design margin on oxidizer load & ullage (20 %)
    thrust_shape="rectangular (F_avg for t_burn)",  # [A] real curve in Fig.5 is nearly flat
    Cd_list=[0.3, 0.5, 0.7],  # [A] drag coefficient sweep (not in paper)
)
R = {"paper[P]": P, "assumptions[A]": A}

# ======== Step 2: consistency checks of the paper's own numbers ========
chk = {}
chk["F*t [N s]"] = P["F_avg"] * P["t_burn"]
chk["It_text [N s]"] = P["It_text"]
m_prop = P["It_text"] / (P["Isp"] * G0)           # [D] assumes Isp = It/(m_prop g0)
m_ox = m_prop - P["m_fuel"]
chk["m_prop = It/(Isp g0) [kg]"] = m_prop
chk["m_ox = m_prop - 64 g [kg]"] = m_ox
chk["O/F implied"] = m_ox / P["m_fuel"]
chk["Fig.2 sucrose O/F range read off"] = "about 1.5 - 3.0"
chk["Table6 wet-dry [kg] (team A)"] = P["A_wet"] - P["A_dry"]
chk["Table7 wet-dry [kg] (team B)"] = P["B_wet"] - P["B_dry"]
mdot = m_prop / P["t_burn"]; mdot_ox = m_ox / P["t_burn"]; mdot_f = P["m_fuel"] / P["t_burn"]
chk["mdot_total [kg/s]"] = mdot
chk["mdot_ox [kg/s]"] = mdot_ox
chk["mdot_fuel [kg/s]"] = mdot_f
# oxidizer mass flux: paper gives 70 (Table2) and 300 (text); compute ours for several areas
for name, d in [("port D=50mm (Table1)", 0.050), ("chamber D=65mm (Fig.1b read-off)", 0.065)]:
    chk[f"G_ox over {name} [kg/m2 s]"] = mdot_ox / (math.pi * d**2 / 4)
chk["port D that would give G=70 [mm]"] = 1e3 * math.sqrt(4 * mdot_ox / (math.pi * 70))
chk["port D that would give G=300 [mm]"] = 1e3 * math.sqrt(4 * mdot_ox / (math.pi * 300))
R["step2_consistency"] = chk

# ======== Step 4: motor sizing ========
def cf_ideal(g, pe_over_pc, pa_over_pc):
    """Ideal thrust coefficient (standard isentropic formula), given pe/pc and pa/pc."""
    t1 = 2 * g**2 / (g - 1) * (2 / (g + 1)) ** ((g + 1) / (g - 1))
    t2 = 1 - pe_over_pc ** ((g - 1) / g)
    return math.sqrt(t1 * t2) + (pe_over_pc - pa_over_pc) * eps_of(g, pe_over_pc)

def eps_of(g, pr):  # expansion ratio Ae/At for given pe/pc
    num = ((2 / (g + 1)) ** (1 / (g - 1))) * (pr ** (1 / g))
    den = math.sqrt((g + 1) / (g - 1) * (1 - pr ** ((g - 1) / g)))
    return 1 / (num * den) * 1  # Ae/At = 1/( ... )

# optimal expansion at sea level: pe = pa
g = A["gamma"]; pa_pc = PA / P["Pc"]
CF_opt = cf_ideal(g, pa_pc, pa_pc)
eps_opt = eps_of(g, pa_pc)
CF_conv = math.sqrt(2 * g**2 / (g - 1) * (2 / (g + 1)) ** ((g + 1) / (g - 1)) * (1 - pa_pc ** ((g - 1) / g)))  # = CF_opt
At = P["F_avg"] / (CF_opt * P["Pc"])
cstar = P["Pc"] * At / mdot
Isp_chk = CF_opt * cstar / G0
sz = dict(
    note="CF is IDEAL (assumed gamma, optimum expansion at sea level); real CF is lower, so real throat is somewhat larger.",
    CF_ideal=CF_opt, eps_opt=eps_opt, At_mm2=At * 1e6, d_throat_mm=1e3 * math.sqrt(4 * At / math.pi),
    d_exit_mm=1e3 * math.sqrt(4 * At * eps_opt / math.pi),
    cstar_implied_m_s=cstar, cstar_note="depends on assumed gamma via CF (circular, not independent)",
    Isp_from_CF_cstar_s=Isp_chk,
    Isp_paper_s=P["Isp"],
    note2="Isp_from_CF_cstar equals paper Isp by construction (At was sized from F, Pc, mdot): it is NOT an independent validation.",
)
# oxidizer tank
m_ox_load = m_ox * A["tank_margin"]
sz["N2O_load_with_margin_kg"] = m_ox_load
sz["tank_volume_L_liquid"] = m_ox_load / A["rho_N2O_liq"] * 1e3
sz["tank_volume_L_with_ullage"] = sz["tank_volume_L_liquid"] * A["tank_margin"]
# fuel grain volume in chamber (for reference): annulus not specified in the paper -> only chamber bore volume
sz["chamber_bore_volume_L (82 outer, 50 inner, 339 long, [P])"] = math.pi * 0.050**2 / 4 * 0.339 * 1e3
sz["T/W of candy motor alone not meaningful; see flight T/W below"] = None
R["step4_sizing"] = sz

# ======== Step 5: 1-D flight ========
def isa_rho(h):  # ISA troposphere (standard formula)
    T = 288.15 - 0.0065 * h
    p = 101325 * (T / 288.15) ** 5.25588
    return p / (287.053 * T)

def fly(m_dry_total, m_prop_total, D, Cd, F=P["F_avg"], tb=P["t_burn"], drag=True, dt=1e-3):
    S = math.pi * D**2 / 4
    t = 0.0; h = 0.0; v = 0.0; m = m_dry_total + m_prop_total
    mdot_l = m_prop_total / tb
    tt, hh, vv = [], [], []
    while True:
        burning = t < tb
        thrust = F if burning else 0.0
        D_ = 0.5 * isa_rho(h) * v * abs(v) * Cd * S if drag else 0.0
        # on the pad thrust must exceed weight, otherwise it does not lift off
        a = (thrust - D_) / m - G0
        if h <= 0 and v <= 0 and a < 0: a = 0.0; v = 0.0
        v += a * dt; h += v * dt; t += dt
        if burning: m -= mdot_l * dt
        tt.append(t); hh.append(h); vv.append(v)
        if v < 0 and not burning: break
        if t > 120: break
    i = int(np.argmax(hh))
    return dict(apogee_m=hh[i], t_apogee_s=tt[i], v_burnout_m_s=vv[int(tb / dt)], t=tt, h=hh, v=vv)

def summarize(name, m_dry, mp, D, alt_paper):
    """Thrust is scaled with propellant mass (same Isp & burn time): F = Isp*g0*mp/tb."""
    F = P["Isp"] * G0 * mp / P["t_burn"]
    m_prop_ = mp
    out = {"m_dry": m_dry, "m_prop": mp, "D": D, "paper_calc_alt_m": alt_paper, "F_used_N": F}
    r0 = fly(m_dry, mp, D, 0, F=F, drag=False)
    out["drag_free_apogee_m"] = r0["apogee_m"]; out["drag_free_t_apogee_s"] = r0["t_apogee_s"]
    out["TW_liftoff"] = F / ((m_dry + mp) * G0)
    for cd in A["Cd_list"]:
        r = fly(m_dry, mp, D, cd, F=F)
        out[f"apogee_Cd{cd}_m"] = r["apogee_m"]; out[f"t_apogee_Cd{cd}_s"] = r["t_apogee_s"]
    # inferred Cd that would reproduce the paper's calculated altitude (INFERENCE, only if reachable)
    lo, hi = 0.0, 5.0
    if r0["apogee_m"] >= alt_paper:
        for _ in range(40):
            mid = (lo + hi) / 2
            if fly(m_dry, mp, D, mid, F=F)["apogee_m"] > alt_paper: lo = mid
            else: hi = mid
        out["Cd_needed_to_match_paper_calc"] = (lo + hi) / 2
    else:
        out["Cd_needed_to_match_paper_calc"] = "UNREACHABLE: even drag-free apogee is below the paper's number"
    return out, r0

fl = {}
# Candy rocket (Table 5): dry 10 kg. Whether "dry" includes the empty motor is not stated -> two readings.
cand, r_candy_nodrag = summarize("candy", P["rocket_dry"], m_prop, P["rocket_D"], P["alt_calc_candy"])
fl["candy_rocket_dry10kg+prop"] = cand
cand2, _ = summarize("candy", P["rocket_dry"] - m_prop, m_prop, P["rocket_D"], P["alt_calc_candy"])
fl["candy_rocket_10kg_is_launch_mass(wet)"] = cand2
# Team A: wet-dry = 570 g as propellant
a, _ = summarize("A", P["A_dry"], P["A_wet"] - P["A_dry"], P["A_D"], P["A_alt_calc"])
fl["teamA_(Table6)"] = a
# Team B: wet-dry = 330 g -> only about 58 % of candy-motor propellant. Flag.
b, _ = summarize("B", P["B_dry"], P["B_wet"] - P["B_dry"], P["B_D"], P["B_alt_calc"])
fl["teamB_(Table7)"] = b
# What impulse would team B need to reach 1187 m drag-free? (inverse problem)
def drag_free_apogee_for_prop(m_dry, mp, tb=P["t_burn"], F_scale=None):
    F = (P["Isp"] * G0 * mp) / tb
    return fly(m_dry, mp, 0.089, 0, F=F, tb=tb, drag=False)["apogee_m"]
need = None
for mp in np.arange(0.2, 1.5, 0.005):
    if drag_free_apogee_for_prop(P["B_dry"], mp) >= P["B_alt_calc"]:
        need = float(mp); break
fl["teamB_propellant_needed_for_1187m_drag_free_kg_(Isp168.2, same burn time)"] = need
# inverse: impulse needed to reach the paper's apogees drag-free with 10 kg DRY (face value of Table 5)
inv = {}
for tgt in (313.0, 319.0):
    for mp in np.arange(0.3, 1.5, 0.001):
        F_ = P["Isp"] * G0 * mp / P["t_burn"]
        if fly(P["rocket_dry"], mp, P["rocket_D"], 0, F=F_, drag=False)["apogee_m"] >= tgt:
            inv[f"target_{int(tgt)}m"] = dict(m_prop_kg=float(mp), It_Ns=float(P["Isp"] * G0 * mp)); break
fl["inverse_dragfree_10kg_dry"] = inv
# stoichiometric O/F for sucrose + N2O -> CO2 + H2O + N2  (C12H22O11 + 24 N2O -> 12 CO2 + 11 H2O + 24 N2)
M_suc = 12*12.011 + 22*1.008 + 11*15.999; M_n2o = 2*14.007 + 15.999
chk["stoichiometric O/F (sucrose+N2O)"] = 24 * M_n2o / M_suc
chk["implied O/F / stoich"] = chk["O/F implied"] / chk["stoichiometric O/F (sucrose+N2O)"]
R["step5_flight"] = fl

# ---------------- plots ----------------
plt.rcParams["font.family"] = ["WenQuanYi Zen Hei", "DejaVu Sans"]
fig, ax = plt.subplots(1, 2, figsize=(11, 4))
t = np.array(r_candy_nodrag["t"]); h = np.array(r_candy_nodrag["h"]); v = np.array(r_candy_nodrag["v"])
ax[0].plot(t, h, label="空気抵抗なし (drag-free)")
for cd in A["Cd_list"]:
    r = fly(P["rocket_dry"], m_prop, P["rocket_D"], cd)
    ax[0].plot(r["t"], r["h"], label=f"Cd={cd} [仮定]")
ax[0].axhline(P["alt_meas_candy"], ls="--", c="k", lw=0.8); ax[0].text(0.3, P["alt_meas_candy"] + 5, "論文の実測 319 m", fontsize=8); ax[0].set_xlim(0, 12)
ax[0].set_xlabel("時間 [s]"); ax[0].set_ylabel("高度 [m]"); ax[0].legend(fontsize=8); ax[0].set_title("キャンディロケット(Table 5)の飛行計算")
ax[1].bar(["図2の上限(3)", "実験の逆算"], [3.0, chk["O/F implied"]], color=["tab:blue", "tab:red"])
ax[1].set_ylabel("O/F (酸化剤質量 ÷ 燃料質量)"); ax[1].set_title("混合比のズレ")
plt.tight_layout(); plt.savefig(OUT / "flight_and_OF.png", dpi=150)

(OUT / "results.json").write_text(json.dumps(R, indent=2, ensure_ascii=False, default=str))
print(json.dumps({k: R[k] for k in ["step2_consistency", "step4_sizing", "step5_flight"]},
                 indent=2, ensure_ascii=False, default=str))
