#!/usr/bin/env python3
"""Convert a (steady, pre-converged) case to a transient continuation with
buoyantBoussinesqPimpleFoam, restarting from its latest steady iteration.

Buoyant room flow is intrinsically unsteady, so the steady SIMPLE run is used
only as an initial condition. The transient run spins up for SPIN seconds and
then time-averages for AVG seconds (~1.2 nominal air-change periods, V/Q = 411 s)."""
import os, re, sys

SPIN, AVG = float(os.environ.get("SPIN", 240)), float(os.environ.get("AVG", 480))


def convert(case):
    ts = [float(d) for d in os.listdir(case) if re.fullmatch(r"[0-9.]+", d)]
    t0 = max(ts)
    t0s = "%g" % t0
    # steady averages must not be carried over into the transient average
    for f in ("UMean", "TMean", "ageMean", "kMean"):
        p = os.path.join(case, t0s, f)
        if os.path.exists(p):
            os.remove(p)
    cd = open(os.path.join(case, "system/controlDict")).read()
    head, funcs = cd.split("functions", 1)
    tend = t0 + SPIN + AVG
    head = f"""FoamFile {{ version 2.0; format ascii; class dictionary; object controlDict; }}
application buoyantBoussinesqPimpleFoam;
startFrom latestTime; startTime {t0s}; stopAt endTime; endTime {tend:g};
deltaT 0.05; adjustTimeStep yes; maxCo 5; maxDeltaT 0.5;
writeControl adjustable; writeInterval {tend - t0:g}; purgeWrite 0;
writeFormat ascii; writePrecision 7; writeCompression off;
timeFormat general; timePrecision 8; runTimeModifiable true;

"""
    funcs = re.sub(r"timeStart \S+;", f"timeStart {t0 + SPIN:g}; restartOnRestart true;", funcs)
    funcs = funcs.replace("writeControl timeStep; writeInterval 20;",
                          "writeControl adjustable; writeInterval 5;")
    open(os.path.join(case, "system/controlDict"), "w").write(head + "functions" + funcs)
    fs = open(os.path.join(case, "system/fvSchemes")).read()
    fs = fs.replace("default steadyState;", "default backward;")
    fs = fs.replace("div(phi,U) bounded Gauss", "div(phi,U) Gauss")
    fs = fs.replace("div(phi,T) bounded Gauss", "div(phi,T) Gauss")
    fs = fs.replace("div(phi,age) bounded Gauss", "div(phi,age) Gauss")
    fs = fs.replace("div(phi,k) bounded Gauss", "div(phi,k) Gauss")
    fs = fs.replace("div(phi,epsilon) bounded Gauss", "div(phi,epsilon) Gauss")
    open(os.path.join(case, "system/fvSchemes"), "w").write(fs)
    open(os.path.join(case, "system/fvSolution"), "w").write("""FoamFile { version 2.0; format ascii; class dictionary; object fvSolution; }
solvers
{
    p_rgh { solver GAMG; smoother DIC; tolerance 1e-7; relTol 0.01; }
    p_rghFinal { $p_rgh; relTol 0; }
    "(U|T|k|epsilon)" { solver PBiCGStab; preconditioner DILU; tolerance 1e-8; relTol 0.01; }
    "(U|T|k|epsilon)Final" { $U; relTol 0; }
    age { solver PBiCGStab; preconditioner DILU; tolerance 1e-8; relTol 0; }
}
PIMPLE
{
    momentumPredictor yes;
    nOuterCorrectors 1;
    nCorrectors 2;
    nNonOrthogonalCorrectors 0;
    pRefCell 0; pRefValue 0;
}
relaxationFactors { equations { ".*" 1; } }
""")
    open(os.path.join(case, "Allrun.transient"), "w").write("""#!/bin/bash
cd "${0%/*}" || exit 1
. /usr/lib/openfoam/openfoam2506/etc/bashrc
rm -f done
buoyantBoussinesqPimpleFoam > log.transient 2>&1 &&
postProcess -func writeCellCentres -latestTime > log.cc 2>&1 &&
postProcess -func writeCellVolumes -latestTime > log.cv 2>&1
echo "$? $(basename $PWD)" > done
""")
    os.chmod(os.path.join(case, "Allrun.transient"), 0o755)


if __name__ == "__main__":
    for c in sys.argv[1:]:
        convert(c)
        print("converted", c)
