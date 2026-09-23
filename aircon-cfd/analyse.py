#!/usr/bin/env python3
"""Post-process the living-room cases: comfort / efficiency metrics + figures."""
import glob, json, math, os, re, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import make_cases as mc

HERE = os.path.dirname(os.path.abspath(__file__))
CASES = sys.argv[1] if len(sys.argv) > 1 else mc.CASES
OUT = os.environ.get("OUT", os.path.join(HERE, "results"))
os.makedirs(OUT, exist_ok=True)

V_ROOM = mc.LX * mc.LY * mc.LZ
TAU_N = V_ROOM / mc.Q                     # nominal time constant, s
SET = {"cool": 24.0, "heat": 21.0}        # occupied-zone setpoints

INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#898781", "#e1e0d9"
CMAP = LinearSegmentedColormap.from_list(
    "div", ["#104281", "#2a78d6", "#86b6ef", "#f0efec", "#f3b08f", "#eb6834", "#a23a12"])
plt.rcParams.update({"font.size": 9, "axes.edgecolor": "#c3c2b7", "axes.labelcolor": INK2,
                     "xtick.color": MUTED, "ytick.color": MUTED, "axes.titlecolor": INK,
                     "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb"})


def read_field(path):
    s = open(path).read()
    m = re.search(r"internalField\s+nonuniform\s+List<(\w+)>\s*(\d+)\s*\(", s)
    n, i = int(m.group(2)), m.end()
    body = s[i:s.index("\n)", i)]
    if m.group(1) == "vector":
        a = np.array(body.replace("(", " ").replace(")", " ").split(), float)
        return a.reshape(n, 3)
    return np.array(body.split(), float)


def last_time(case):
    ts = [d for d in os.listdir(case) if re.fullmatch(r"[0-9.]+", d) and d != "0"]
    return max(ts, key=float)


def avg_start(case):
    s = open(os.path.join(case, "system/controlDict")).read()
    return float(re.search(r"timeStart (\S+);", s).group(1))


def fo_mean(case, fo, col):
    """Mean/std of a function-object signal over the fieldAverage window
    (the transient run; the steady spin-up writes to a different time dir)."""
    fs = sorted(glob.glob(os.path.join(case, "postProcessing", fo, "*", "*.dat")),
                key=lambda p: float(p.split(os.sep)[-2]))
    a = np.atleast_2d(np.loadtxt(fs[-1], comments="#"))
    w = a[a[:, 0] >= avg_start(case)]
    return w[:, col].mean(), w[:, col].std(), a


def meta(case):
    return dict(l.strip().split("=") for l in open(os.path.join(case, "meta.txt")))


def cop_rel(mode, Tsup):
    """Carnot-fraction COP (constant fraction cancels in ratios)."""
    if mode == "cool":
        Te, Tc = Tsup - 7.0 + 273.15, 35.0 + 12.0 + 273.15   # 35 C outdoor
        return Te / (Tc - Te)
    Tc, Te = Tsup + 8.0 + 273.15, 7.0 - 7.0 + 273.15         # 7 C outdoor
    return Tc / (Tc - Te)


def analyse(case):
    md = meta(case)
    mode, t = md["mode"], last_time(case)
    ld = lambda f: read_field(os.path.join(case, t, f))
    C, V = ld("C"), ld("V")
    T, U, age = ld("TMean"), ld("UMean"), ld("ageMean")
    x, y, z = C.T
    spd = np.linalg.norm(U, axis=1)
    (lo, hi) = mc.ZONES["occ"]
    occ = (x > lo[0]) & (x < hi[0]) & (y > lo[1]) & (y < hi[1]) & (z > lo[2]) & (z < hi[2])
    Vo = V[occ]
    Tocc = np.average(T[occ], weights=Vo)
    Tsup = float(md["Tsup"])
    Tret, Tret_sd, _ = fo_mean(case, "Tret", 1)
    mcp = mc.RHOCP * mc.Q                                   # W/K of supply air
    ua = "UA" in md
    if ua:   # heating: envelope losses UA_i (T_zone_i - Tout)
        Tout = float(md["Tout"])
        def zoneT(zn):
            lo_, hi_ = mc.ZONES[zn]
            sel = (x > lo_[0]) & (x < hi_[0]) & (y > lo_[1]) & (y < hi_[1]) & (z > lo_[2]) & (z < hi_[2])
            return np.average(T[sel], weights=V[sel])
        loss = sum(u * (zoneT(zn) - Tout) for zn, u in mc.MODES["heat"]["loads"].items())
        Tret_bal = Tsup - loss / mcp
    else:
        Tret_bal = Tsup + float(md["load"]) / mcp
    # temperature effectiveness from simultaneous window averages of the
    # return and occupied zone (self-consistent even if a small global drift
    # remains; Tret_bal - Tret is reported as the energy-balance residual)
    eps = (Tret - Tsup) / (Tocc - Tsup)
    # vertical profile over the occupied footprint
    foot = (x > lo[0]) & (x < hi[0]) & (y > lo[1]) & (y < hi[1])
    zs = np.unique(np.round(z[foot], 4))
    prof = np.array([np.average(T[foot & (np.abs(z - zz) < 1e-3)],
                                weights=V[foot & (np.abs(z - zz) < 1e-3)]) for zz in zs])
    Tz = lambda h: float(np.interp(h, zs, prof))
    # ADPI (EDT in K and m/s form)
    edt = (T[occ] - Tocc) - 7.66 * (spd[occ] - 0.15)
    ok = (edt > -1.7) & (edt < 1.1) & (spd[occ] < 0.35)
    adpi = 100 * Vo[ok].sum() / Vo.sum()
    draught = 100 * Vo[spd[occ] > 0.25].sum() / Vo.sum()
    age_occ = np.average(age[occ], weights=Vo)
    age_ret, _, _ = fo_mean(case, "Tret", 2)
    # energy to hold the occupied zone at the setpoint (room sensor), same airflow
    sp = SET[mode]
    if ua:
        # no fixed sources -> with the flow pattern frozen, (T - Tout) scales
        # linearly with (Tsup - Tout); this keeps the extra envelope loss of a
        # hot ceiling layer in the bill
        kk = (Tocc - Tout) / (Tsup - Tout)
        Tsup_req = Tout + (sp - Tout) / kk
        heat_req = mcp * (Tsup - Tret) * (Tsup_req - Tout) / (Tsup - Tout)
        UAt = float(md["UA"])
        heat_ideal = UAt * (sp - Tout)
        Tsup_ideal = sp + heat_ideal / mcp
        elec_vs_ideal = (heat_req / cop_rel(mode, Tsup_req)) / (heat_ideal / cop_rel(mode, Tsup_ideal))
        dTu = heat_req / mcp
    else:
        dTu = float(md["load"]) / mcp
        Tsup_req = sp - dTu / eps
        Tsup_ideal = sp - dTu
        heat_req = heat_ideal = float(md["load"])
        elec_vs_ideal = cop_rel(mode, Tsup_ideal) / cop_rel(mode, Tsup_req)
    r = dict(case=os.path.basename(case), geom=md["geom"], mode=mode, angle=int(md["angle"]),
             refine=int(md["refine"]), vsup=float(md["vsup"]),
             Tsup=Tsup, Tret_fo=Tret, Tret_fo_sd=Tret_sd, Tret_bal=Tret_bal,
             Tocc=Tocc, eps=eps, bypass_equiv=1 - eps,
             ret_minus_occ=Tret - Tocc, balance_residual=Tret_bal - Tret,
             T_ankle=Tz(0.1), T_seated_head=Tz(1.1), T_stand_head=Tz(1.7), T_ceiling=Tz(2.8),
             strat_01_11=Tz(1.1) - Tz(0.1),
             spd_occ=np.average(spd[occ], weights=Vo), draught_pct=draught, adpi=adpi,
             age_occ=age_occ, age_ret=age_ret, tau_n=TAU_N, ace=age_ret / age_occ,
             Tsup_req=Tsup_req, heat_req=heat_req, heat_ideal=heat_ideal,
             elec_vs_ideal=elec_vs_ideal)
    return r, dict(C=C, T=T, U=U, spd=spd, age=age, zs=zs, prof=prof, md=md)


def grid(C):
    xs, ys, zs = (np.unique(np.round(C[:, i], 5)) for i in range(3))
    return xs, ys, zs


def slice_plots(name, d, r):
    C, T, U = d["C"], d["T"], d["U"]
    xs, ys, zs = grid(C)
    nx, ny, nz = len(xs), len(ys), len(zs)
    T3, U3 = T.reshape(nz, ny, nx), U.reshape(nz, ny, nx, 3)
    geom = mc.GEOMS[d["md"]["geom"]]
    sx = 0.5 * sum(geom["sup"])
    mode = d["md"]["mode"]
    vmin, vmax = (16, 28) if mode == "cool" else (17, 29)
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.3), gridspec_kw=dict(width_ratios=[4.6, 4.1, 4.1]))
    # (a) section through supply centre (y-z)
    i = np.argmin(abs(xs - sx))
    ax = axs[0]
    cf = ax.contourf(ys, zs, T3[:, :, i], levels=np.linspace(vmin, vmax, 25), cmap=CMAP, extend="both")
    sk = (slice(None, None, 3), slice(None, None, 2))
    ax.quiver(ys[::2], zs[::3], U3[:, :, i, 1][sk], U3[:, :, i, 2][sk], color=INK, alpha=.55,
              scale=12, width=.0022)
    ax.set_title(f"Section through supply (x = {xs[i]:.2f} m)", loc="left")
    ax.set_xlabel("distance from unit wall  y [m]"); ax.set_ylabel("height z [m]")
    ax.axhspan(0.1, 1.8, xmin=0.3/4.6, xmax=3.6/4.6, fill=False, ls="--", ec=INK2, lw=.8)
    # (b) section across room at mid depth (x-z)
    j = np.argmin(abs(ys - 2.3))
    ax = axs[1]
    ax.contourf(xs, zs, T3[:, j, :], levels=np.linspace(vmin, vmax, 25), cmap=CMAP, extend="both")
    ax.quiver(xs[::2], zs[::3], U3[:, j, :, 0][sk], U3[:, j, :, 2][sk], color=INK, alpha=.55,
              scale=12, width=.0022)
    ax.set_title(f"Cross-section at mid-room (y = {ys[j]:.2f} m)", loc="left")
    ax.set_xlabel("x [m]")
    # (c) plan at grille height, near the unit wall (short-circuit zone)
    k = np.argmin(abs(zs - 2.69))
    ax = axs[2]
    ax.contourf(xs, ys, T3[k], levels=np.linspace(vmin, vmax, 25), cmap=CMAP, extend="both")
    sk2 = (slice(None, None, 2), slice(None, None, 2))
    ax.quiver(xs[::2], ys[::2], U3[k, :, :, 0][sk2], U3[k, :, :, 1][sk2], color=INK, alpha=.55,
              scale=12, width=.0022)
    ax.plot(geom["sup"], [0.03, 0.03], color="#2a78d6" if mode == "cool" else "#eb6834", lw=5,
            solid_capstyle="butt", label="supply")
    rr = geom["ret"]
    if rr[0] == "wall":
        ax.plot(rr[1:3], [0.03, 0.03], color=INK2, lw=5, solid_capstyle="butt", label="return")
    else:
        ax.add_patch(plt.Rectangle((rr[1], rr[3]), rr[2] - rr[1], rr[4] - rr[3], fill=False,
                                   ec=INK2, lw=2, label="return (ceiling)"))
    ax.legend(loc="upper right", fontsize=8, frameon=False)
    ax.set_title(f"Plan at grille height (z = {zs[k]:.2f} m)", loc="left")
    ax.set_xlabel("x [m]  (along unit wall)"); ax.set_ylabel("y [m]")
    for a in axs:
        a.set_aspect("equal")
    cb = fig.colorbar(cf, ax=axs, shrink=.85, pad=.01)
    cb.set_label("mean air temperature [°C]")
    fig.suptitle(f"{name}   ·   occupied-zone {r['Tocc']:.1f} °C, return {r['Tret_fo']:.1f} °C, "
                 f"ankle→head ΔT {r['strat_01_11']:+.1f} K, ε = {r['eps']:.2f}",
                 x=0.01, ha="left", color=INK, fontsize=10.5)
    fig.savefig(os.path.join(OUT, f"slices_{name}.png"), dpi=110, bbox_inches="tight")
    plt.close(fig)


LABEL = {"A1_adjacent": "A1 twin grilles, adjacent",
         "A2_separated": "A2 twin grilles, opposite ends",
         "A3_fullwidth": "A3 full-width split strip",
         "B_ceiling": "B wall supply + ceiling return"}
COL = {"A1_adjacent": "#eb6834", "A2_separated": "#1baf7a",
       "A3_fullwidth": "#eda100", "B_ceiling": "#2a78d6"}


def profile_plot(data):
    fig, axs = plt.subplots(1, 2, figsize=(10, 4.6), sharey=True)
    for ax, mode in zip(axs, ("cool", "heat")):
        for name, (r, d) in sorted(data.items()):
            if r["mode"] != mode or r["refine"] != 1:
                continue
            ls = "--" if r["angle"] else "-"
            ax.plot(d["prof"], d["zs"], ls, color=COL[r["geom"]], lw=2,
                    label=LABEL[r["geom"]] + (f", vanes {r['angle']}° down" if r["angle"] else ""))
        ax.axhspan(0.1, 1.8, color=GRID, alpha=.5, lw=0)
        ax.grid(color=GRID, lw=.6)
        ax.set_title("Cooling (supply 13 °C)" if mode == "cool" else "Heating (supply 34 °C)", loc="left")
        ax.set_xlabel("mean temperature over occupied footprint [°C]")
    axs[0].set_ylabel("height [m]  (shaded = occupied zone)")
    axs[1].legend(fontsize=7.5, frameon=False, loc="lower right")
    axs[0].legend(fontsize=7.5, frameon=False, loc="lower right")
    fig.savefig(os.path.join(OUT, "profiles.png"), dpi=120, bbox_inches="tight")
    plt.close(fig)


def main():
    data = {}
    for case in sorted(glob.glob(os.path.join(CASES, "*"))):
        if not os.path.exists(os.path.join(case, "done")) or open(os.path.join(case, "done")).read()[0] != "0":
            print("skip", case); continue
        r, d = analyse(case)
        data[r["case"]] = (r, d)
        slice_plots(r["case"], d, r)
        print(f"{r['case']:28s} Tocc {r['Tocc']:5.2f} Tret {r['Tret_bal']:5.2f}/{r['Tret_fo']:5.2f}±{r['Tret_fo_sd']:.2f} "
              f"eps {r['eps']:.3f} strat {r['strat_01_11']:+.2f} ADPI {r['adpi']:4.0f} "
              f"draught {r['draught_pct']:4.0f}% ACE {r['ace']:.2f} elec {r['elec_vs_ideal']:.3f}")
    profile_plot(data)
    json.dump([r for r, _ in data.values()], open(os.path.join(OUT, "metrics.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
