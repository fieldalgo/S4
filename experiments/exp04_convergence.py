"""What a converged answer looks like, and how to tell you have not got one.

R+T must be 1 for a lossless structure, so the departure from 1 is a free
convergence indicator that needs no reference.  The requested and the actual
order count differ under the default circular truncation, which is worth seeing
once: asking for 100 gets 97, and asking for a rectangular set gets all of them.
"""
import S4

NAME = "convergence, and requested vs actual orders"


def solve(nb, truncation=None):
    S = S4.New(Lattice=((1, 0), (0, 1)), NumBasis=nb)
    opts = dict(PolarizationDecomposition=True, PolarizationBasis="Normal")
    if truncation:
        opts["LatticeTruncation"] = truncation
    S.SetOptions(**opts)
    S.SetMaterial(Name="bg", Epsilon=1.0)
    S.SetMaterial(Name="si", Epsilon=3.5 ** 2)
    S.AddLayer(Name="in", Thickness=0, Material="bg")
    S.AddLayer(Name="grid", Thickness=0.25, Material="bg")
    S.SetRegionCircle(Layer="grid", Material="si", Center=(0, 0), Radius=0.3)
    S.AddLayer(Name="out", Thickness=0, Material="si")
    S.SetExcitationPlanewave(IncidenceAngles=(0, 0), sAmplitude=1, pAmplitude=0)
    S.SetFrequency(1 / 1.55)
    inc, back = S.GetPowerFlux(Layer="in")
    fwd, _ = S.GetPowerFlux(Layer="out")
    i = complex(inc).real
    o = S.GetOptions()
    return o, complex(fwd).real / i, abs(complex(back).real) / i


def run():
    rows = []
    for nb in (25, 100, 200, 400, 800):
        o, T, R = solve(nb)
        rows.append(["%d" % nb, "%d" % o["NumBasis"], "%.9f" % T,
                     "%+.1e" % (T + R - 1)])
    last = float(rows[-1][3])
    trunc = []
    for t in ("Circular", "Parallelogramic"):
        o, T, _ = solve(289, t)
        trunc.append([t, "%d" % o["NumBasisRequested"], "%d" % o["NumBasis"],
                      "%.9f" % T])
    return {
        "headers": ["NumBasis asked", "actual orders", "T", "R+T-1"],
        "rows": rows,
        "extra": {
            "headers": ["truncation", "requested", "actual", "T"],
            "rows": trunc,
            "caption": "The same request, two truncations. GetOptions reports "
                       "both counts; the difference is how a truncation that "
                       "was not applied shows itself.",
        },
        "summary": "R+T-1 falls from %s to %s over NumBasis 25 to 800"
                   % (rows[0][3], rows[-1][3]),
        "pass": abs(last) < 1e-3,
    }
