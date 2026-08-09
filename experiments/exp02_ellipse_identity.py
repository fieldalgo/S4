"""An ellipse with equal half-widths is a circle, at every angle.

Two identities that need no reference implementation: an ellipse whose
half-widths are equal is a circle and rotating it changes nothing, and an
ellipse rotated ninety degrees is the same shape with its half-widths swapped.
The circle reaches the solver through a different branch of the same function,
so the first is an independent route to the same number.
"""
import S4

NAME = "ellipse(r,r) is a circle"
CENTRE = (0.05, -0.03)          # off-centre, so a rotation has something to do


def transmission(place, basis, nb=100):
    S = S4.New(Lattice=((1, 0), (0, 1)), NumBasis=nb)
    S.SetOptions(PolarizationDecomposition=True, PolarizationBasis=basis)
    S.SetMaterial(Name="bg", Epsilon=1.0)
    S.SetMaterial(Name="rod", Epsilon=3.5 ** 2)
    S.AddLayer(Name="top", Thickness=0, Material="bg")
    S.AddLayer(Name="pat", Thickness=0.15, Material="bg")
    place(S)
    S.AddLayer(Name="bot", Thickness=0, Material="bg")
    S.SetExcitationPlanewave(IncidenceAngles=(20, 30), sAmplitude=1, pAmplitude=0)
    S.SetFrequency(1 / 1.6)
    inc, _ = S.GetPowerFlux(Layer="top")
    fwd, _ = S.GetPowerFlux(Layer="bot")
    return complex(fwd).real / complex(inc).real


def run():
    rows, worst = [], 0.0
    for basis in ("Normal", "Jones", "Default"):
        circle = transmission(
            lambda S: S.SetRegionCircle(Layer="pat", Material="rod",
                                        Center=CENTRE, Radius=0.30), basis)
        for angle in (0, 25, 45, 90, 143):
            got = transmission(
                lambda S, a=angle: S.SetRegionEllipse(
                    Layer="pat", Material="rod", Center=CENTRE, Angle=a,
                    Halfwidths=(0.30, 0.30)), basis)
            d = abs(got - circle)
            worst = max(worst, d)
            rows.append([basis, "%d" % angle, "%.12f" % circle,
                         "%.12f" % got, "%.1e" % d])
    return {
        "headers": ["basis", "angle", "circle T", "ellipse(r,r) T", "error"],
        "rows": rows,
        "summary": "worst departure from the circle %.1e across three bases "
                   "and five angles" % worst,
        "pass": worst < 1e-11,
    }
