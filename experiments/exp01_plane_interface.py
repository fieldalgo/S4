"""A plane interface, against the closed form.

Two layers and nothing between them is the simplest structure S4 can express
and the only one with an analytic answer.  It returned T = 0, R = 0, R+T = 0
before the two-layer fix; the point of measuring it here is that Fresnel is an
answer from outside S4, so agreement is evidence rather than self-consistency.
"""
import cmath

import S4

NAME = "plane interface vs Fresnel"


def fresnel(n2, theta_deg, pol):
    th = cmath.pi * theta_deg / 180
    c1, s1 = cmath.cos(th), cmath.sin(th)
    c2 = cmath.sqrt(1 - (s1 / n2) ** 2)
    r = ((c1 - n2 * c2) / (c1 + n2 * c2) if pol == "s"
         else (n2 * c1 - c2) / (n2 * c1 + c2))
    return 1 - abs(r) ** 2, abs(r) ** 2


def s4_interface(n2, theta, pol):
    S = S4.New(Lattice=((1, 0), (0, 1)), NumBasis=9)
    S.SetMaterial(Name="a", Epsilon=1.0)
    S.SetMaterial(Name="b", Epsilon=complex(n2) ** 2)
    S.AddLayer(Name="top", Thickness=0, Material="a")
    S.AddLayer(Name="bot", Thickness=0, Material="b")
    s_amp, p_amp = (1, 0) if pol == "s" else (0, 1)
    S.SetExcitationPlanewave(IncidenceAngles=(theta, 0),
                             sAmplitude=s_amp, pAmplitude=p_amp)
    S.SetFrequency(1 / 1.55)
    inc, back = S.GetPowerFlux(Layer="top")
    fwd, _ = S.GetPowerFlux(Layer="bot")
    i = complex(inc).real
    return complex(fwd).real / i, abs(complex(back).real) / i


def run():
    rows, worst = [], 0.0
    for n2 in (1.5, 3.5):
        for theta in (0, 20, 45, 70):
            for pol in ("s", "p"):
                T, R = s4_interface(n2, theta, pol)
                wT, wR = fresnel(n2, theta, pol)
                d = max(abs(T - wT.real), abs(R - wR.real))
                worst = max(worst, d)
                rows.append((n2, theta, pol, T, wT.real, d, T + R))
    return {
        "headers": ["n2", "angle", "pol", "S4 T", "Fresnel T", "error", "R+T"],
        "rows": [["%.1f" % r[0], "%d" % r[1], r[2], "%.9f" % r[3],
                  "%.9f" % r[4], "%.1e" % r[5], "%.9f" % r[6]] for r in rows],
        "summary": "worst disagreement with Fresnel %.1e over %d cases"
                   % (worst, len(rows)),
        "pass": worst < 1e-12,
    }
