#!/usr/bin/env python3
"""Generate OpenFOAM (v1912, buoyantBoussinesqSimpleFoam) cases for the
living-room ducted air-conditioner study.

Room: 4.1 m (x, unit wall width) x 4.6 m (y, depth) x 2.9 m (z, height).
Unit wall is y = 0. Grilles: top 150 mm below ceiling, 120 mm tall
(z = 2.63 .. 2.75). The supply grille is modelled as a 60 mm effective slot
(z = 2.66 .. 2.72) - i.e. a ~50% free-area grille - so that both the mass flow
and the jet momentum are right.
"""
import math, os, shutil, sys

HERE = os.path.dirname(os.path.abspath(__file__))
CASES = os.path.join(HERE, "cases")

LX, LY, LZ = 4.1, 4.6, 2.9
Q = 0.133                    # m3/s, Daikin FDXM35 nominal airflow
RHOCP = 1.2 * 1005.0         # J/m3K
Z_GRILLE = (2.63, 2.75)      # visual grille (return)
Z_SUPPLY = (2.66, 2.72)      # effective supply slot (50% free area)

# ---------------------------------------------------------------- geometries
# x-extents are on the 0.1 m mesh lines
GEOMS = {
    # symmetrical twin grilles side by side, 100 mm apart
    "A1_adjacent":  dict(sup=(0.8, 2.0), ret=("wall", 2.1, 3.3)),
    # symmetrical twin grilles pushed to opposite ends (1.3 m gap)
    "A2_separated": dict(sup=(0.2, 1.4), ret=("wall", 2.7, 3.9)),
    # one continuous full-width strip: left half supply, right half return
    "A3_fullwidth": dict(sup=(0.1, 2.0), ret=("wall", 2.1, 4.0)),
    # reference fix for heating: same high supply, return low on the same wall
    "C_lowreturn":  dict(sup=(0.8, 2.0), ret=("lowwall", 2.1, 3.3)),
    # centred wall supply + 600x400 ceiling return near the far wall
    "B_ceiling":    dict(sup=(1.4, 2.6), ret=("ceiling", 1.7, 2.3, 3.7, 4.1)),
}

# ------------------------------------------------------------ thermal modes
# Heat sources (W) injected into thin cell layers (exact energy input,
# independent of wall functions). Radiation is not modelled.
MODES = {
    "cool": dict(Tsup=13.0, Tinit=24.0, TRef=24.0, loads={
        "farWall": 1000.0,   # glazed / sun-exposed far wall
        "floor":    400.0,   # solar patch on floor + furniture
        "ceiling":  150.0,   # roof gain
        "people":   250.0,   # 2 people + TV etc. in occupied zone
    }),
    # Heating: envelope losses are UA*(T_local - T_out) (semi-implicit), so a
    # cold pool loses less than a warm ceiling layer - fixed-watt losses would
    # grossly exaggerate stratification. UA sized for ~1.7 kW at 21 C / 7 C.
    "heat": dict(Tsup=34.0, Tinit=22.0, TRef=22.0, Tout=7.0, loads={
        "farWall": 57.0,     # W/K, window wall
        "sideL":   14.0,
        "sideR":   14.0,
        "floor":   18.0,
        "ceiling": 18.0,
    }),
}

ZONES = {  # boxes for cellZones
    "farWall": ((0, LY - 0.1, 0), (LX, LY, LZ)),
    "sideL":   ((0, 0, 0), (0.1, LY, LZ)),
    "sideR":   ((LX - 0.1, 0, 0), (LX, LY, LZ)),
    "floor":   ((0, 0, 0), (LX, LY, 0.1)),
    "ceiling": ((0, 0, LZ - 0.03), (LX, LY, LZ)),
    "people":  ((1.0, 1.5, 0.0), (3.0, 3.5, 1.2)),
    # ASHRAE 55 occupied zone: 0.1-1.8 m high, 1 m off the window wall
    "occ":     ((0.3, 0.3, 0.1), (3.8, 3.6, 1.8)),
}

# ------------------------------------------------------------------- cases
# name: (geometry, mode, supply vane angle below horizontal, mesh refine)
RUNS = {
    "cool_A1_adjacent":      ("A1_adjacent",  "cool", 0, 1),
    "cool_A2_separated":     ("A2_separated", "cool", 0, 1),
    "cool_A3_fullwidth":     ("A3_fullwidth", "cool", 0, 1),
    "cool_B_ceiling":        ("B_ceiling",    "cool", 0, 1),
    "heat_A1_adjacent":      ("A1_adjacent",  "heat", 0, 1),
    "heat_A2_separated":     ("A2_separated", "heat", 0, 1),
    "heat_A3_fullwidth":     ("A3_fullwidth", "heat", 0, 1),
    "heat_B_ceiling":        ("B_ceiling",    "heat", 0, 1),
    "heat_A1_adjacent_down30":  ("A1_adjacent",  "heat", 30, 1),
    "heat_A2_separated_down30": ("A2_separated", "heat", 30, 1),
    "heat_B_ceiling_down30":    ("B_ceiling",    "heat", 30, 1),
    "cool_A1_adjacent_down30":  ("A1_adjacent",  "cool", 30, 1),
    "heat_C_lowreturn":         ("C_lowreturn",  "heat", 0, 1),
    "heat_A1_adjacent_down45":  ("A1_adjacent",  "heat", 45, 1),
    "heat_A2_separated_down45": ("A2_separated", "heat", 45, 1),
    # mesh-sensitivity checks (2x finer grille / ceiling-jet band)
    "cool_A1_adjacent_fine": ("A1_adjacent",  "cool", 0, 2),
    "heat_A1_adjacent_fine": ("A1_adjacent",  "heat", 0, 2),
}

HDR = """FoamFile
{
    version     2.0;
    format      ascii;
    class       %s;
    object      %s;
}
"""


def w(case, rel, cls, body):
    p = os.path.join(case, rel)
    os.makedirs(os.path.dirname(p), exist_ok=True)
    with open(p, "w") as f:
        f.write(HDR % (cls, os.path.basename(rel)) + body)


def vec(v):
    return "(%g %g %g)" % tuple(v)


def field(case, name, cls, dim, internal, bcs):
    body = "dimensions %s;\ninternalField %s;\nboundaryField\n{\n" % (dim, internal)
    for patch, bc in bcs.items():
        body += "    %s\n    {\n%s    }\n" % (patch, "".join(
            "        %s %s;\n" % kv for kv in bc.items()))
    body += "}\n"
    w(case, "0/" + name, cls, body)


def make(name, geom, mode, angle, refine):
    g, m = GEOMS[geom], MODES[mode]
    case = os.path.join(CASES, name)
    shutil.rmtree(case, ignore_errors=True)
    os.makedirs(case)

    # ---------------- mesh: dx = dy = 0.1; dz = 0.1 below 2.3 m and
    # 0.03/refine m in the top 0.6 m (grille / ceiling-jet band)
    nx, ny, nzt = 41, 46, 20 * refine
    w(case, "system/blockMeshDict", "dictionary", f"""
convertToMeters 1;
vertices
(
 (0 0 0) ({LX} 0 0) ({LX} {LY} 0) (0 {LY} 0)
 (0 0 {LZ}) ({LX} 0 {LZ}) ({LX} {LY} {LZ}) (0 {LY} {LZ})
);
blocks
(
 hex (0 1 2 3 4 5 6 7) ({nx} {ny} {23 + nzt})
 simpleGrading (1 1 ((2.3 23 1) (0.6 {nzt} 1)))
);
boundary
(
 walls {{ type wall; faces ((0 4 7 3) (1 2 6 5) (0 1 5 4) (3 7 6 2) (0 3 2 1) (4 5 6 7)); }}
);
""")

    # ---------------- patches + zones
    e = 1e-3
    sx0, sx1 = g["sup"]
    actions = [f"""{{ name supplyF; type faceSet; action new; source boxToFace;
   sourceInfo {{ box ({sx0-e} {-e} {Z_SUPPLY[0]-e}) ({sx1+e} {e} {Z_SUPPLY[1]+e}); }} }}"""]
    r = g["ret"]
    if r[0] == "wall":
        rbox = f"({r[1]-e} {-e} {Z_GRILLE[0]-e}) ({r[2]+e} {e} {Z_GRILLE[1]+e})"
        ret_area = (r[2] - r[1]) * (Z_GRILLE[1] - Z_GRILLE[0])
    elif r[0] == "lowwall":   # z = 0.1 .. 0.2 (on the 0.1 m mesh lines)
        rbox = f"({r[1]-e} {-e} {0.1-e}) ({r[2]+e} {e} {0.2+e})"
        ret_area = (r[2] - r[1]) * 0.1
    else:
        rbox = f"({r[1]-e} {r[3]-e} {LZ-e}) ({r[2]+e} {r[4]+e} {LZ+e})"
        ret_area = (r[2] - r[1]) * (r[4] - r[3])
    actions.append(f"""{{ name returnF; type faceSet; action new; source boxToFace;
   sourceInfo {{ box {rbox}; }} }}""")
    for z, (lo, hi) in ZONES.items():
        actions.append(f"""{{ name {z}; type cellSet; action new; source boxToCell;
   sourceInfo {{ box {vec(lo)} {vec(hi)}; }} }}""")
        actions.append(f"""{{ name {z}Zone; type cellZoneSet; action new; source setToCellZone;
   sourceInfo {{ set {z}; }} }}""")
    w(case, "system/topoSetDict", "dictionary",
      "actions\n(\n" + "\n".join(actions) + "\n);\n")
    w(case, "system/createPatchDict", "dictionary", """
pointSync false;
patches
(
 { name supply; patchInfo { type patch; } constructFrom set; set supplyF; }
 { name return; patchInfo { type patch; } constructFrom set; set returnF; }
);
""")

    # ---------------- supply velocity (normal flux exact, vanes tilt it down)
    sup_area = (sx1 - sx0) * (Z_SUPPLY[1] - Z_SUPPLY[0])
    vn = Q / sup_area
    Usup = (0.0, vn, -vn * math.tan(math.radians(angle)))

    # ---------------- constant
    w(case, "constant/g", "uniformDimensionedVectorField",
      "dimensions [0 1 -2 0 0 0 0];\nvalue (0 0 -9.81);\n")
    w(case, "constant/transportProperties", "dictionary", f"""
transportModel Newtonian;
nu 1.5e-05;
beta 3.4e-03;
TRef {m['TRef']};
Pr 0.71;
Prt 0.85;
""")
    w(case, "constant/turbulenceProperties", "dictionary", """
simulationType RAS;
RAS { RASModel RNGkEpsilon; turbulence on; printCoeffs on; }
""")
    fvo = ""
    for zone, W in m["loads"].items():
        if "Tout" in m:   # UA [W/K]: S = UA*Tout - UA*T
            su, sp = W * m["Tout"] / RHOCP, -W / RHOCP
        else:             # fixed heat gain [W]
            su, sp = W / RHOCP, 0.0
        fvo += f"""
{zone}Load
{{
    type scalarSemiImplicitSource;
    active yes;
    scalarSemiImplicitSourceCoeffs
    {{
        selectionMode cellZone;
        cellZone {zone}Zone;
        volumeMode absolute;
        injectionRateSuSp {{ T ({su:.6g} {sp:.6g}); }}
    }}
}}
"""
    w(case, "constant/fvOptions", "dictionary", fvo)

    # ---------------- 0/
    Ts, T0 = m["Tsup"], m["Tinit"]
    wall = "walls"
    field(case, "U", "volVectorField", "[0 1 -1 0 0 0 0]", "uniform (0 0 0)", {
        wall: {"type": "noSlip"},
        "supply": {"type": "fixedValue", "value": "uniform " + vec(Usup)},
        "return": {"type": "pressureInletOutletVelocity", "value": "uniform (0 0 0)"},
    })
    field(case, "p_rgh", "volScalarField", "[0 2 -2 0 0 0 0]", "uniform 0", {
        wall: {"type": "fixedFluxPressure", "rho": "rhok", "value": "uniform 0"},
        "supply": {"type": "fixedFluxPressure", "rho": "rhok", "value": "uniform 0"},
        "return": {"type": "fixedValue", "value": "uniform 0"},
    })
    field(case, "p", "volScalarField", "[0 2 -2 0 0 0 0]", "uniform 0", {
        p: {"type": "calculated", "value": "uniform 0"} for p in (wall, "supply", "return")})
    field(case, "T", "volScalarField", "[0 0 0 1 0 0 0]", f"uniform {T0}", {
        wall: {"type": "zeroGradient"},
        "supply": {"type": "fixedValue", "value": f"uniform {Ts}"},
        "return": {"type": "inletOutlet", "inletValue": f"uniform {T0}", "value": f"uniform {T0}"},
    })
    field(case, "age", "volScalarField", "[0 0 0 0 0 0 0]", "uniform 0", {
        wall: {"type": "zeroGradient"},
        "supply": {"type": "fixedValue", "value": "uniform 0"},
        "return": {"type": "inletOutlet", "inletValue": "uniform 0", "value": "uniform 0"},
    })
    field(case, "k", "volScalarField", "[0 2 -2 0 0 0 0]", "uniform 1e-3", {
        wall: {"type": "kqRWallFunction", "value": "uniform 1e-3"},
        "supply": {"type": "turbulentIntensityKineticEnergyInlet", "intensity": "0.05", "value": "uniform 1e-3"},
        "return": {"type": "inletOutlet", "inletValue": "uniform 1e-3", "value": "uniform 1e-3"},
    })
    field(case, "epsilon", "volScalarField", "[0 2 -3 0 0 0 0]", "uniform 1e-4", {
        wall: {"type": "epsilonWallFunction", "value": "uniform 1e-4"},
        "supply": {"type": "turbulentMixingLengthDissipationRateInlet", "mixingLength": "0.008", "value": "uniform 1e-4"},
        "return": {"type": "inletOutlet", "inletValue": "uniform 1e-4", "value": "uniform 1e-4"},
    })
    field(case, "nut", "volScalarField", "[0 2 -1 0 0 0 0]", "uniform 0", {
        wall: {"type": "nutkWallFunction", "value": "uniform 0"},
        "supply": {"type": "calculated", "value": "uniform 0"},
        "return": {"type": "calculated", "value": "uniform 0"},
    })
    field(case, "alphat", "volScalarField", "[0 2 -1 0 0 0 0]", "uniform 0", {
        wall: {"type": "alphatJayatillekeWallFunction", "Prt": "0.85", "value": "uniform 0"},
        "supply": {"type": "calculated", "value": "uniform 0"},
        "return": {"type": "calculated", "value": "uniform 0"},
    })

    # ---------------- system
    iters, avg_from = 1000, 500
    w(case, "system/controlDict", "dictionary", f"""
application buoyantBoussinesqSimpleFoam;
startFrom latestTime; startTime 0; stopAt endTime; endTime {iters};
deltaT 1; writeControl timeStep; writeInterval {iters}; purgeWrite 0;
writeFormat ascii; writePrecision 7; writeCompression off;
timeFormat general; timePrecision 6; runTimeModifiable true;

functions
{{
    age
    {{
        type scalarTransport; libs (solverFunctionObjects);
        field age; bounded01 false; alphaD 1; alphaDt 1;
        writeControl writeTime;
        fvOptions
        {{
            unitySource
            {{
                type scalarSemiImplicitSource; enabled true;
                scalarSemiImplicitSourceCoeffs
                {{
                    selectionMode all; volumeMode specific;
                    injectionRateSuSp {{ age (1 0); }}
                }}
            }}
        }}
    }}
    avg
    {{
        type fieldAverage; libs (fieldFunctionObjects);
        timeStart {avg_from}; writeControl writeTime;
        fields ( U {{ mean on; prime2Mean off; base iteration; }}
                 T {{ mean on; prime2Mean off; base iteration; }}
                 age {{ mean on; prime2Mean off; base iteration; }}
                 k {{ mean on; prime2Mean off; base iteration; }} );
    }}
    Tret
    {{
        type surfaceFieldValue; libs (fieldFunctionObjects);
        writeControl timeStep; writeInterval 20; log false; writeFields false;
        regionType patch; name return; operation weightedAverage; weightField phi;
        fields (T age);
    }}
    Tocc
    {{
        type volFieldValue; libs (fieldFunctionObjects);
        writeControl timeStep; writeInterval 20; log false; writeFields false;
        regionType cellZone; name occZone; operation volAverage;
        fields (T age);
    }}
}}
""")
    w(case, "system/fvSchemes", "dictionary", """
ddtSchemes { default steadyState; }
gradSchemes { default Gauss linear; grad(U) cellLimited Gauss linear 1; }
divSchemes
{
    default none;
    div(phi,U) bounded Gauss linearUpwindV grad(U);
    div(phi,T) bounded Gauss linearUpwind grad(T);
    div(phi,age) bounded Gauss upwind;
    div(phi,k) bounded Gauss upwind;
    div(phi,epsilon) bounded Gauss upwind;
    div((nuEff*dev2(T(grad(U))))) Gauss linear;
}
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
wallDist { method meshWave; }
""")
    w(case, "system/fvSolution", "dictionary", """
solvers
{
    p_rgh { solver GAMG; smoother GaussSeidel; tolerance 1e-7; relTol 0.01; }
    "(U|T|k|epsilon|age)" { solver PBiCGStab; preconditioner DILU; tolerance 1e-7; relTol 0.1; }
}
SIMPLE
{
    nNonOrthogonalCorrectors 0;
    pRefCell 0; pRefValue 0;
}
relaxationFactors
{
    fields { p_rgh 0.5; }
    equations { U 0.3; T 0.5; age 0.7; "(k|epsilon)" 0.5; }
}
""")
    with open(os.path.join(case, "meta.txt"), "w") as f:
        f.write(f"geom={geom}\nmode={mode}\nangle={angle}\nrefine={refine}\n"
                f"sup_area={sup_area}\nret_area={ret_area}\nvsup={vn}\n"
                f"Tsup={Ts}\nQ={Q}\nload={sum(m['loads'].values())}\n"
                + (f"Tout={m['Tout']}\nUA={sum(m['loads'].values())}\n" if "Tout" in m else ""))
    with open(os.path.join(case, "Allrun"), "w") as f:
        f.write("""#!/bin/bash
cd "${0%/*}" || exit 1
. /usr/lib/openfoam/openfoam2506/etc/bashrc
blockMesh > log.blockMesh 2>&1 &&
topoSet > log.topoSet 2>&1 &&
createPatch -overwrite > log.createPatch 2>&1 &&
checkMesh > log.checkMesh 2>&1 &&
buoyantBoussinesqSimpleFoam > log.solver 2>&1 &&
python3 ../../to_transient.py . > log.convert 2>&1 &&
exec ./Allrun.transient
echo "$? $(basename $PWD)" > done
""")
    os.chmod(os.path.join(case, "Allrun"), 0o755)


if __name__ == "__main__":
    sel = sys.argv[1:] or list(RUNS)
    for n in sel:
        make(n, *RUNS[n])
        print("made", n)
