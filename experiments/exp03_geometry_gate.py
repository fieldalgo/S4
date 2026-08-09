"""Does S4 hold the structure it was handed?

The mean permittivity over the cell is the DC Fourier coefficient, and the
reconstruction is band-limited, so sampling it past Nyquist gives that
coefficient exactly -- no truncation error, at any basis size.  It therefore
separates "the geometry went in wrong" from "it has not converged yet", which
is the distinction most of the defects in this fork hid behind.

Worth running first on any new structure: an error here is not a convergence
problem and no amount of NumBasis will move it.
"""
import math

import S4

NAME = "geometry gate: mean permittivity"
EPS_BG, EPS_ROD = 1.0, 3.5 ** 2


def mean_eps(place, n=64, thick=0.2):
    S = S4.New(Lattice=((1, 0), (0, 1)), NumBasis=25)
    S.SetMaterial(Name="bg", Epsilon=EPS_BG)
    S.SetMaterial(Name="rod", Epsilon=EPS_ROD)
    S.AddLayer(Name="top", Thickness=0, Material="bg")
    S.AddLayer(Name="pat", Thickness=thick, Material="bg")
    place(S)
    S.AddLayer(Name="bot", Thickness=0, Material="bg")
    S.SetExcitationPlanewave(IncidenceAngles=(0, 0), sAmplitude=1, pAmplitude=0)
    S.SetFrequency(1 / 1.55)
    S.GetPowerFlux(Layer="top")
    z = thick / 2
    return sum(S.GetEpsilon((i + 0.5) / n - 0.5, (j + 0.5) / n - 0.5, z).real
               for i in range(n) for j in range(n)) / (n * n)


def run():
    d = EPS_ROD - EPS_BG
    cases = [
        ("circle r=0.30",
         lambda S: S.SetRegionCircle(Layer="pat", Material="rod",
                                     Center=(0, 0), Radius=0.30),
         EPS_BG + d * math.pi * 0.30 ** 2),
        ("square 0.5 x 0.5",
         lambda S: S.SetRegionRectangle(Layer="pat", Material="rod",
                                        Center=(0, 0), Angle=0,
                                        Halfwidths=(0.25, 0.25)),
         EPS_BG + d * 0.25),
        ("ellipse 0.34 x 0.18 at 25 deg",
         lambda S: S.SetRegionEllipse(Layer="pat", Material="rod",
                                      Center=(0, 0), Angle=25,
                                      Halfwidths=(0.34, 0.18)),
         EPS_BG + d * math.pi * 0.34 * 0.18),
        ("circle in a ring (nested)",
         lambda S: (S.SetRegionCircle(Layer="pat", Material="rod",
                                      Center=(0, 0), Radius=0.30),
                    S.SetRegionCircle(Layer="pat", Material="bg",
                                      Center=(0, 0), Radius=0.15)),
         EPS_BG + d * math.pi * (0.30 ** 2 - 0.15 ** 2)),
    ]
    rows, worst = [], 0.0
    for label, place, exact in cases:
        got = mean_eps(place)
        err = abs(got - exact)
        worst = max(worst, err)
        rows.append([label, "%.9f" % exact, "%.9f" % got, "%.1e" % err])
    return {
        "headers": ["structure", "exact", "S4", "error"],
        "rows": rows,
        "summary": "worst geometric error %.1e; this is exact at any NumBasis"
                   % worst,
        "pass": worst < 1e-9,
    }
