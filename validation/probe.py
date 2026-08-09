"""Exercise every entry point of the S4 Python binding and report numbers.

Each test returns a flat dict of name -> real or complex number.  Nothing in
here asserts: the point is to emit values that can be compared digit for digit
against a different build of the same library, and to make every binding
function actually execute so that a port which merely links is not mistaken for
a port that works.

Run one test and print its results as JSON:

    python probe.py run <test-name>

List the tests:

    python probe.py list

The driver in compare.py runs each test in its own process, so a crash in one
test is recorded rather than taking the whole run down with it.
"""

import json
import math
import sys
import traceback

import s4compat
from s4compat import New, Sim

TESTS = {}


def test(fn):
    TESTS[fn.__name__] = fn
    return fn


# --------------------------------------------------------------------------
# shared geometry: the benchmark this project has used all along -- a square
# lattice of dielectric rods, period 1, radius 0.3, n = 3.5, thickness 0.15,
# lambda 1.6, normal incidence, s-polarised.
# --------------------------------------------------------------------------

def rod_sim(nb=200, basis="Normal", radius=0.30, thick=0.15):
    S = New(Lattice=((1, 0), (0, 1)), NumBasis=nb)
    if basis is not None:
        S.SetOptions(PolarizationDecomposition=True, PolarizationBasis=basis)
    S.SetMaterial(Name="air", Epsilon=1.0)
    S.SetMaterial(Name="rod", Epsilon=3.5 ** 2)
    S.AddLayer(Name="top", Thickness=0, Material="air")
    S.AddLayer(Name="pat", Thickness=thick, Material="air")
    S.SetRegionCircle(Layer="pat", Material="rod", Center=(0, 0), Radius=radius)
    S.AddLayer(Name="bot", Thickness=0, Material="air")
    S.SetExcitationPlanewave(IncidenceAngles=(0, 0), sAmplitude=1, pAmplitude=0)
    S.SetFrequency(1 / 1.6)
    return S


def flux(S, layer="bot", top="top"):
    inc, back = S.GetPowerFlux(Layer=top)
    fwd, _ = S.GetPowerFlux(Layer=layer)
    return fwd.real / inc.real, abs(back.real) / inc.real


# --------------------------------------------------------------------------
# 1. the established reference numbers
# --------------------------------------------------------------------------

@test
def reference_convergence():
    """The four NumBasis values whose T is recorded to nine digits."""
    out = {}
    for nb in (200, 400, 600, 800):
        T, R = flux(rod_sim(nb))
        out["T_nb%d" % nb] = T
        out["RpT_nb%d" % nb] = R + T
    return out


@test
def basis_variants():
    """The four Fourier-factorisation settings at NumBasis 200."""
    out = {}
    for tag, basis in (("none", None), ("default", "Default"),
                       ("normal", "Normal"), ("jones", "Jones")):
        T, R = flux(rod_sim(200, basis))
        out["T_" + tag] = T
        out["RpT_" + tag] = R + T
    return out


# --------------------------------------------------------------------------
# 2. all four region types
# --------------------------------------------------------------------------

def _patterned(nb, place, basis="Normal"):
    S = New(Lattice=((1, 0), (0, 1)), NumBasis=nb)
    if basis is not None:
        S.SetOptions(PolarizationDecomposition=True, PolarizationBasis=basis)
    S.SetMaterial(Name="air", Epsilon=1.0)
    S.SetMaterial(Name="rod", Epsilon=3.5 ** 2)
    S.AddLayer(Name="top", Thickness=0, Material="air")
    S.AddLayer(Name="pat", Thickness=0.15, Material="air")
    place(S)
    S.AddLayer(Name="bot", Thickness=0, Material="air")
    S.SetExcitationPlanewave(IncidenceAngles=(0, 0), sAmplitude=1, pAmplitude=0)
    S.SetFrequency(1 / 1.6)
    return S


@test
def region_circle_vs_ellipse():
    """A circle and a same-radius ellipse must give bit-comparable answers."""
    c = _patterned(120, lambda S: S.SetRegionCircle(
        Layer="pat", Material="rod", Center=(0, 0), Radius=0.30))
    e = _patterned(120, lambda S: S.SetRegionEllipse(
        Layer="pat", Material="rod", Center=(0, 0), Angle=0,
        Halfwidths=(0.30, 0.30)))
    Tc, Rc = flux(c)
    Te, Re = flux(e)
    return {"T_circle": Tc, "T_ellipse": Te, "diff": abs(Tc - Te),
            "RpT_circle": Rc + Tc, "RpT_ellipse": Re + Te}


@test
def region_rectangle():
    """A rectangle, and the same rectangle expressed as a polygon."""
    r = _patterned(120, lambda S: S.SetRegionRectangle(
        Layer="pat", Material="rod", Center=(0, 0), Angle=0,
        Halfwidths=(0.25, 0.15)))
    p = _patterned(120, lambda S: S.SetRegionPolygon(
        Layer="pat", Material="rod", Center=(0, 0), Angle=0,
        Vertices=((-0.25, -0.15), (0.25, -0.15), (0.25, 0.15), (-0.25, 0.15))))
    Tr, Rr = flux(r)
    Tp, Rp = flux(p)
    return {"T_rect": Tr, "T_poly": Tp, "diff": abs(Tr - Tp),
            "RpT_rect": Rr + Tr, "RpT_poly": Rp + Tp}


@test
def region_rectangle_rotated():
    """Rotating a rectangle by 90 degrees must equal swapping its halfwidths.

    The rectangle is centred off the origin, so this also drives the rotation
    branch of shape_contains_point with a non-trivial centre.
    """
    a = _patterned(120, lambda S: S.SetRegionRectangle(
        Layer="pat", Material="rod", Center=(0.11, 0.23), Angle=90,
        Halfwidths=(0.25, 0.15)))
    b = _patterned(120, lambda S: S.SetRegionRectangle(
        Layer="pat", Material="rod", Center=(0.11, 0.23), Angle=0,
        Halfwidths=(0.15, 0.25)))
    Ta, Ra = flux(a)
    Tb, Rb = flux(b)
    return {"T_rot90": Ta, "T_swapped": Tb, "diff": abs(Ta - Tb),
            "RpT_rot90": Ra + Ta, "RpT_swapped": Rb + Tb}


@test
def region_polygon_offcentre():
    """A polygon at a general centre and angle -- pure execution + energy."""
    verts = ((-0.20, -0.10), (0.22, -0.14), (0.18, 0.19), (-0.12, 0.21))
    S = _patterned(120, lambda S: S.SetRegionPolygon(
        Layer="pat", Material="rod", Center=(0.07, -0.05), Angle=33,
        Vertices=verts))
    T, R = flux(S)
    return {"T": T, "RpT": R + T}


@test
def region_multiple_disjoint():
    """Four disjoint rods -- exercises the containment tree with siblings."""
    def place(S):
        for cx, cy in ((-0.25, -0.25), (0.25, -0.25), (-0.25, 0.25), (0.25, 0.25)):
            S.SetRegionCircle(Layer="pat", Material="rod",
                              Center=(cx, cy), Radius=0.12)
    S = _patterned(120, place)
    T, R = flux(S)
    return {"T": T, "RpT": R + T}


@test
def region_nested():
    """A ring: a big rod with a small air rod inside it (proper nesting)."""
    def place(S):
        S.SetRegionCircle(Layer="pat", Material="rod",
                          Center=(0, 0), Radius=0.35)
        S.SetRegionCircle(Layer="pat", Material="air",
                          Center=(0, 0), Radius=0.15)
    S = _patterned(120, place)
    T, R = flux(S)
    return {"T": T, "RpT": R + T}


@test
def region_removal():
    """RemoveLayerRegions must restore the unpatterned layer exactly."""
    S = _patterned(80, lambda S: S.SetRegionCircle(
        Layer="pat", Material="rod", Center=(0, 0), Radius=0.30))
    S.RemoveLayerRegions(Layer="pat")
    T, R = flux(S)

    U = New(Lattice=((1, 0), (0, 1)), NumBasis=80)
    U.SetOptions(PolarizationDecomposition=True, PolarizationBasis="Normal")
    U.SetMaterial(Name="air", Epsilon=1.0)
    U.AddLayer(Name="top", Thickness=0, Material="air")
    U.AddLayer(Name="pat", Thickness=0.15, Material="air")
    U.AddLayer(Name="bot", Thickness=0, Material="air")
    U.SetExcitationPlanewave(IncidenceAngles=(0, 0), sAmplitude=1, pAmplitude=0)
    U.SetFrequency(1 / 1.6)
    Tu, Ru = flux(U)
    return {"T_removed": T, "T_never_patterned": Tu, "diff": abs(T - Tu),
            "RpT_removed": R + T}


# --------------------------------------------------------------------------
# 3. layer copy -- the path whose guard was broken by the integer-handle port
# --------------------------------------------------------------------------

@test
def layer_copy():
    """A copied patterned layer must behave as a second identical layer."""
    def build(use_copy):
        S = New(Lattice=((1, 0), (0, 1)), NumBasis=120)
        S.SetOptions(PolarizationDecomposition=True, PolarizationBasis="Normal")
        S.SetMaterial(Name="air", Epsilon=1.0)
        S.SetMaterial(Name="rod", Epsilon=3.5 ** 2)
        S.AddLayer(Name="top", Thickness=0, Material="air")
        S.AddLayer(Name="pat", Thickness=0.15, Material="air")
        S.SetRegionCircle(Layer="pat", Material="rod",
                          Center=(0, 0), Radius=0.30)
        S.AddLayer(Name="gap", Thickness=0.22, Material="air")
        if use_copy:
            S.AddLayerCopy(Name="pat2", Thickness=0.15, Layer="pat")
        else:
            S.AddLayer(Name="pat2", Thickness=0.15, Material="air")
            S.SetRegionCircle(Layer="pat2", Material="rod",
                              Center=(0, 0), Radius=0.30)
        S.AddLayer(Name="bot", Thickness=0, Material="air")
        S.SetExcitationPlanewave(IncidenceAngles=(0, 0), sAmplitude=1, pAmplitude=0)
        S.SetFrequency(1 / 1.6)
        return S

    Tc, Rc = flux(build(True))
    Te, Re = flux(build(False))
    return {"T_copy": Tc, "T_explicit": Te, "diff": abs(Tc - Te),
            "RpT_copy": Rc + Tc}


@test
def layer_copy_guard():
    """Patterning must be refused on a copy and allowed on an ordinary layer.

    Before b6fe76b the guard read a S4_LayerID as if it were a pointer, so it
    fired for every ordinary layer (id 0 and up) and never for a copy (id -1).
    The two booleans below are what that bug inverted.
    """
    S = New(Lattice=((1, 0), (0, 1)), NumBasis=20)
    S.SetMaterial(Name="air", Epsilon=1.0)
    S.SetMaterial(Name="rod", Epsilon=3.5 ** 2)
    S.AddLayer(Name="pat", Thickness=0.15, Material="air")
    S.AddLayerCopy(Name="cp", Thickness=0.15, Layer="pat")

    ordinary_ok = 1.0
    try:
        S.SetRegionCircle(Layer="pat", Material="rod",
                          Center=(0, 0), Radius=0.2)
    except Exception:
        ordinary_ok = 0.0

    copy_refused = 0.0
    try:
        S.SetRegionCircle(Layer="cp", Material="rod",
                          Center=(0, 0), Radius=0.2)
    except Exception:
        copy_refused = 1.0

    # the patterning above must have actually landed: a patterned layer has a
    # position-dependent epsilon, an unpatterned one does not.
    eps_in = S.GetEpsilon(0.0, 0.0, 0.05)
    eps_out = S.GetEpsilon(0.45, 0.45, 0.05)
    return {"ordinary_layer_patternable": ordinary_ok,
            "copy_layer_refused": copy_refused,
            "eps_at_centre": eps_in, "eps_at_corner": eps_out}


@test
def set_layer_reuse():
    """SetLayer on an existing name must redefine it; on a new name, create it."""
    S = New(Lattice=((1, 0), (0, 1)), NumBasis=80)
    S.SetMaterial(Name="air", Epsilon=1.0)
    S.SetMaterial(Name="glass", Epsilon=2.25)
    S.AddLayer(Name="top", Thickness=0, Material="air")
    S.AddLayer(Name="mid", Thickness=0.5, Material="air")
    S.SetLayer(Name="mid", Thickness=0.5, Material="glass")   # redefine
    S.SetLayer(Name="bot", Thickness=0, Material="air")       # create
    S.SetExcitationPlanewave(IncidenceAngles=(0, 0), sAmplitude=1, pAmplitude=0)
    S.SetFrequency(1 / 1.5)
    T, R = flux(S)

    n, d, lam = 1.5, 0.5, 1.5
    k = 2 * math.pi / lam
    r = (1 - n) / (1 + n)
    ph = math.cos(2 * n * k * d), math.sin(2 * n * k * d)
    num_re = r * (1 - ph[0])
    num_im = r * (-ph[1])
    den_re = 1 - r * r * ph[0]
    den_im = r * r * ph[1]
    den2 = den_re * den_re + den_im * den_im
    rr = ((num_re * den_re + num_im * den_im) ** 2 +
          (num_im * den_re - num_re * den_im) ** 2) / (den2 * den2)
    return {"T": T, "RpT": R + T, "R_fresnel": rr, "R_err": abs(R - rr)}


@test
def set_layer_thickness():
    """Changing thickness through SetLayerThickness must equal building it so."""
    def build(t, via_setter):
        S = New(Lattice=((1, 0), (0, 1)), NumBasis=80)
        S.SetMaterial(Name="air", Epsilon=1.0)
        S.SetMaterial(Name="glass", Epsilon=2.25)
        S.AddLayer(Name="top", Thickness=0, Material="air")
        S.AddLayer(Name="mid", Thickness=0.11 if via_setter else t,
                   Material="glass")
        S.AddLayer(Name="bot", Thickness=0, Material="air")
        if via_setter:
            S.SetLayerThickness(Layer="mid", Thickness=t)
        S.SetExcitationPlanewave(IncidenceAngles=(0, 0), sAmplitude=1, pAmplitude=0)
        S.SetFrequency(1 / 1.31)
        return flux(S)[0]

    a = build(0.37, True)
    b = build(0.37, False)
    return {"T_via_setter": a, "T_direct": b, "diff": abs(a - b)}


# --------------------------------------------------------------------------
# 4. diffraction orders
# --------------------------------------------------------------------------

@test
def lattice_truncation():
    """LatticeTruncation must actually choose the basis it names.

    G_select runs inside S4_Simulation_New, before SetOptions can be called, so
    on a build that never re-selects the option is accepted, stored, and has no
    effect: the basis stays circular whatever was asked for.  Nothing reports
    it, and the count can even come out right by coincidence while the set is
    different -- NumBasis 441 and 1369 both do.

    Parallelogramic with NumBasis (2k+1)^2 has one right answer: the full
    (2k+1) x (2k+1) block of orders.  That is what makes NumBasis comparable
    with a code that truncates rectangularly.
    """
    out = {}
    for nb, half in ((289, 8), (441, 10), (625, 12)):
        want = set((a, b) for a in range(-half, half + 1)
                   for b in range(-half, half + 1))
        for tag, trunc in (("par", "Parallelogramic"), ("cir", "Circular")):
            S = New(Lattice=((2.5, 0), (0, 2.5)), NumBasis=nb)
            S.SetOptions(LatticeTruncation=trunc)
            S.SetMaterial(Name="a", Epsilon=1.0)
            S.AddLayer(Name="t", Thickness=0, Material="a")
            S.AddLayer(Name="b", Thickness=0, Material="a")
            S.SetExcitationPlanewave(IncidenceAngles=(0, 0), sAmplitude=1,
                                     pAmplitude=0)
            S.SetFrequency(1 / 1.0)
            S.GetPowerFlux(Layer="t")
            G = set((int(g[0]), int(g[1])) for g in S.GetBasisSet())
            out["n_%s_%d" % (tag, nb)] = float(len(G))
            out["is_rect_%s_%d" % (tag, nb)] = 1.0 if G == want else 0.0
    # the two truncations must not agree; if they do, one of them is ignored
    for nb in (289, 441, 625):
        out["differ_%d" % nb] = (1.0 if out["is_rect_par_%d" % nb] !=
                                 out["is_rect_cir_%d" % nb] else 0.0)
    return out


@test
def basis_set():
    """GetBasisSet must return exactly the G vectors the solve uses."""
    S = rod_sim(120)
    S.GetPowerFlux(Layer="top")
    G = S.GetBasisSet()
    out = {"count": float(len(G))}
    for i in range(6):
        out["G%d_m" % i] = float(G[i][0])
        out["G%d_n" % i] = float(G[i][1])
    # the basis must be symmetric about the origin and contain (0,0) once
    s = set((int(g[0]), int(g[1])) for g in G)
    out["has_zero"] = 1.0 if (0, 0) in s else 0.0
    out["symmetric"] = 1.0 if all((-m, -n) in s for (m, n) in s) else 0.0
    out["distinct"] = float(len(s))
    return out


@test
def orders_sum_to_total():
    """Sum over diffraction orders must reproduce the total power flux."""
    S = rod_sim(120)
    out = {}
    for layer in ("top", "bot"):
        tot_f, tot_b = S.GetPowerFlux(Layer=layer)
        by = S.GetPowerFluxByOrder(Layer=layer)
        sf = sum(o[0].real for o in by)
        sb = sum(o[1].real for o in by)
        out["n_orders_" + layer] = float(len(by))
        out["tot_fwd_" + layer] = tot_f.real
        out["sum_fwd_" + layer] = sf
        out["err_fwd_" + layer] = abs(tot_f.real - sf)
        out["tot_bwd_" + layer] = tot_b.real
        out["sum_bwd_" + layer] = sb
        out["err_bwd_" + layer] = abs(tot_b.real - sb)
    return out


@test
def orders_individual():
    """Per-order transmitted power, for the orders that actually propagate.

    Period 1 at lambda 1.6 leaves only the specular order propagating, so the
    interesting check is a shorter wavelength where several orders open up.
    """
    S = New(Lattice=((1, 0), (0, 1)), NumBasis=120)
    S.SetOptions(PolarizationDecomposition=True, PolarizationBasis="Normal")
    S.SetMaterial(Name="air", Epsilon=1.0)
    S.SetMaterial(Name="rod", Epsilon=3.5 ** 2)
    S.AddLayer(Name="top", Thickness=0, Material="air")
    S.AddLayer(Name="pat", Thickness=0.15, Material="air")
    S.SetRegionCircle(Layer="pat", Material="rod", Center=(0, 0), Radius=0.30)
    S.AddLayer(Name="bot", Thickness=0, Material="air")
    S.SetExcitationPlanewave(IncidenceAngles=(0, 0), sAmplitude=1, pAmplitude=0)
    S.SetFrequency(1 / 0.7)          # lambda 0.7 < period: orders open

    G = S.GetBasisSet()
    inc, _ = S.GetPowerFlux(Layer="top")
    by = S.GetPowerFluxByOrder(Layer="bot")
    out = {"inc": inc.real, "n_orders": float(len(by))}
    want = {(0, 0): "00", (1, 0): "10", (-1, 0): "m10",
            (0, 1): "01", (0, -1): "0m1", (1, 1): "11"}
    for i, g in enumerate(G):
        key = want.get((int(g[0]), int(g[1])))
        if key is not None:
            out["T_" + key] = by[i][0].real / inc.real
    # the four first orders are related by the C4v symmetry of the structure
    out["sym_10_vs_m10"] = abs(out.get("T_10", 0) - out.get("T_m10", 0))
    out["sym_01_vs_0m1"] = abs(out.get("T_01", 0) - out.get("T_0m1", 0))
    tot = sum(o[0].real for o in by) / inc.real
    out["T_total"] = tot
    return out


# --------------------------------------------------------------------------
# 5. oblique and conical incidence
# --------------------------------------------------------------------------

@test
def oblique_incidence():
    """Off-normal incidence in the plane of symmetry, s and p."""
    out = {}
    for theta in (10.0, 30.0, 55.0):
        for pol, (sa, pa) in (("s", (1, 0)), ("p", (0, 1))):
            S = New(Lattice=((1, 0), (0, 1)), NumBasis=120)
            S.SetOptions(PolarizationDecomposition=True,
                         PolarizationBasis="Normal")
            S.SetMaterial(Name="air", Epsilon=1.0)
            S.SetMaterial(Name="rod", Epsilon=3.5 ** 2)
            S.AddLayer(Name="top", Thickness=0, Material="air")
            S.AddLayer(Name="pat", Thickness=0.15, Material="air")
            S.SetRegionCircle(Layer="pat", Material="rod",
                              Center=(0, 0), Radius=0.30)
            S.AddLayer(Name="bot", Thickness=0, Material="air")
            S.SetExcitationPlanewave(IncidenceAngles=(theta, 0),
                                     sAmplitude=sa, pAmplitude=pa)
            S.SetFrequency(1 / 1.6)
            T, R = flux(S)
            tag = "%s_th%g" % (pol, theta)
            out["T_" + tag] = T
            out["RpT_" + tag] = R + T
    return out


@test
def conical_incidence():
    """Non-zero azimuth: the fully conical mounting."""
    out = {}
    for theta, phi in ((30.0, 30.0), (30.0, 45.0), (45.0, 20.0)):
        S = New(Lattice=((1, 0), (0, 1)), NumBasis=120)
        S.SetOptions(PolarizationDecomposition=True, PolarizationBasis="Normal")
        S.SetMaterial(Name="air", Epsilon=1.0)
        S.SetMaterial(Name="rod", Epsilon=3.5 ** 2)
        S.AddLayer(Name="top", Thickness=0, Material="air")
        S.AddLayer(Name="pat", Thickness=0.15, Material="air")
        S.SetRegionCircle(Layer="pat", Material="rod",
                          Center=(0, 0), Radius=0.30)
        S.AddLayer(Name="bot", Thickness=0, Material="air")
        S.SetExcitationPlanewave(IncidenceAngles=(theta, phi),
                                 sAmplitude=0.6, pAmplitude=0.8)
        S.SetFrequency(1 / 1.6)
        T, R = flux(S)
        tag = "th%g_ph%g" % (theta, phi)
        out["T_" + tag] = T
        out["RpT_" + tag] = R + T
    return out


@test
def conical_c4_symmetry():
    """A C4-symmetric structure must give the same T at phi and phi+90."""
    def run(phi):
        S = New(Lattice=((1, 0), (0, 1)), NumBasis=120)
        S.SetOptions(PolarizationDecomposition=True, PolarizationBasis="Normal")
        S.SetMaterial(Name="air", Epsilon=1.0)
        S.SetMaterial(Name="rod", Epsilon=3.5 ** 2)
        S.AddLayer(Name="top", Thickness=0, Material="air")
        S.AddLayer(Name="pat", Thickness=0.15, Material="air")
        S.SetRegionCircle(Layer="pat", Material="rod",
                          Center=(0, 0), Radius=0.30)
        S.AddLayer(Name="bot", Thickness=0, Material="air")
        S.SetExcitationPlanewave(IncidenceAngles=(35.0, phi),
                                 sAmplitude=1, pAmplitude=0)
        S.SetFrequency(1 / 1.6)
        return flux(S)[0]

    a, b = run(15.0), run(105.0)
    return {"T_phi15": a, "T_phi105": b, "diff": abs(a - b)}


# --------------------------------------------------------------------------
# 6. anisotropic / tensor materials
# --------------------------------------------------------------------------

@test
def tensor_diagonal():
    """A diagonal tensor whose entries are equal must equal the scalar case."""
    def build(eps):
        S = New(Lattice=((1, 0), (0, 1)), NumBasis=40)
        S.SetMaterial(Name="air", Epsilon=1.0)
        S.SetMaterial(Name="mat", Epsilon=eps)
        S.AddLayer(Name="top", Thickness=0, Material="air")
        S.AddLayer(Name="mid", Thickness=0.37, Material="mat")
        S.AddLayer(Name="bot", Thickness=0, Material="air")
        S.SetExcitationPlanewave(IncidenceAngles=(0, 0), sAmplitude=1, pAmplitude=0)
        S.SetFrequency(1 / 1.31)
        return flux(S)

    e = 4.0
    Ts, Rs = build(e)
    Tt, Rt = build(((e, 0, 0),
                    (0, e, 0),
                    (0, 0, e)))
    return {"T_scalar": Ts, "T_tensor": Tt, "diff": abs(Ts - Tt),
            "RpT_scalar": Rs + Ts, "RpT_tensor": Rt + Tt}


@test
def tensor_uniaxial():
    """A genuinely anisotropic diagonal tensor: s and p must differ."""
    eps = ((4.00, 0, 0),
           (0, 6.25, 0),
           (0, 0, 4.00))

    def run(sa, pa, theta):
        S = New(Lattice=((1, 0), (0, 1)), NumBasis=40)
        S.SetMaterial(Name="air", Epsilon=1.0)
        S.SetMaterial(Name="mat", Epsilon=eps)
        S.AddLayer(Name="top", Thickness=0, Material="air")
        S.AddLayer(Name="mid", Thickness=0.37, Material="mat")
        S.AddLayer(Name="bot", Thickness=0, Material="air")
        S.SetExcitationPlanewave(IncidenceAngles=(theta, 0),
                                 sAmplitude=sa, pAmplitude=pa)
        S.SetFrequency(1 / 1.31)
        return flux(S)

    Ts, Rs = run(1, 0, 0.0)
    Tp, Rp = run(0, 1, 0.0)
    To, Ro = run(1, 0, 40.0)
    return {"T_s": Ts, "T_p": Tp, "T_s_oblique": To,
            "RpT_s": Rs + Ts, "RpT_p": Rp + Tp, "RpT_s_oblique": Ro + To,
            "s_minus_p": abs(Ts - Tp)}


@test
def tensor_offdiagonal_complex():
    """Off-diagonal complex tensor -- the case that used to SIGSEGV.

    This is the branch that reaches the zgeev workspace query with no rwork
    array of its own.  It is here to be executed, not just to be linked.
    """
    eps = ((complex(4.0, 0.10), complex(0.35, 0.20), 0),
           (complex(-0.35, -0.20), complex(4.0, 0.10), 0),
           (0, 0, complex(4.0, 0.10)))
    S = New(Lattice=((1, 0), (0, 1)), NumBasis=40)
    S.SetMaterial(Name="air", Epsilon=1.0)
    S.SetMaterial(Name="mat", Epsilon=eps)
    S.AddLayer(Name="top", Thickness=0, Material="air")
    S.AddLayer(Name="mid", Thickness=0.37, Material="mat")
    S.AddLayer(Name="bot", Thickness=0, Material="air")
    S.SetExcitationPlanewave(IncidenceAngles=(0, 0), sAmplitude=1, pAmplitude=0)
    S.SetFrequency(1 / 1.31)
    T, R = flux(S)
    return {"T": T, "R": R, "RpT": R + T}


@test
def tensor_offdiagonal_lossless():
    """A Hermitian off-diagonal tensor: lossless, so R+T must be 1.

    eps = [[4, i g, 0], [-i g, 4, 0], [0, 0, 4]] is Hermitian for real g, which
    is the gyrotropic (magneto-optic) form.  Any energy imbalance here is a
    solver defect, not absorption.
    """
    g = 0.4
    eps = ((4.0, complex(0.0, g), 0),
           (complex(0.0, -g), 4.0, 0),
           (0, 0, 4.0))
    out = {}
    for theta in (0.0, 25.0):
        S = New(Lattice=((1, 0), (0, 1)), NumBasis=40)
        S.SetMaterial(Name="air", Epsilon=1.0)
        S.SetMaterial(Name="mat", Epsilon=eps)
        S.AddLayer(Name="top", Thickness=0, Material="air")
        S.AddLayer(Name="mid", Thickness=0.37, Material="mat")
        S.AddLayer(Name="bot", Thickness=0, Material="air")
        S.SetExcitationPlanewave(IncidenceAngles=(theta, 0),
                                 sAmplitude=1, pAmplitude=0)
        S.SetFrequency(1 / 1.31)
        T, R = flux(S)
        out["T_th%g" % theta] = T
        out["RpT_th%g" % theta] = R + T
    return out


@test
def tensor_patterned():
    """An anisotropic inclusion inside a patterned layer."""
    eps = ((9.0, complex(0.5, 0.3), 0),
           (complex(-0.5, -0.3), 12.25, 0),
           (0, 0, 9.0))
    S = New(Lattice=((1, 0), (0, 1)), NumBasis=80)
    S.SetMaterial(Name="air", Epsilon=1.0)
    S.SetMaterial(Name="rod", Epsilon=eps)
    S.AddLayer(Name="top", Thickness=0, Material="air")
    S.AddLayer(Name="pat", Thickness=0.15, Material="air")
    S.SetRegionCircle(Layer="pat", Material="rod", Center=(0, 0), Radius=0.30)
    S.AddLayer(Name="bot", Thickness=0, Material="air")
    S.SetExcitationPlanewave(IncidenceAngles=(0, 0), sAmplitude=1, pAmplitude=0)
    S.SetFrequency(1 / 1.6)
    T, R = flux(S)
    return {"T": T, "R": R, "RpT": R + T}


# --------------------------------------------------------------------------
# 7. field output -- linked but, until now, never called
# --------------------------------------------------------------------------

@test
def fields_uniform_slab():
    """In a uniform structure the field must be the analytic plane wave.

    With no patterning and vacuum everywhere, an s-polarised wave of unit
    amplitude travelling in +z has E = (0, 1, 0) exp(i k z) and H = (-1, 0, 0)
    exp(i k z) in S4's normalisation.  Anything else means the field readout is
    wrong regardless of what the flux says.
    """
    lam = 1.31
    S = New(Lattice=((1, 0), (0, 1)), NumBasis=20)
    S.SetMaterial(Name="air", Epsilon=1.0)
    S.AddLayer(Name="top", Thickness=0, Material="air")
    S.AddLayer(Name="mid", Thickness=0.5, Material="air")
    S.AddLayer(Name="bot", Thickness=0, Material="air")
    S.SetExcitationPlanewave(IncidenceAngles=(0, 0), sAmplitude=1, pAmplitude=0)
    S.SetFrequency(1 / lam)

    out = {}
    for tag, z in (("z0", 0.0), ("z1", 0.17), ("z2", 0.5)):
        E, H = S.GetFields(0.0, 0.0, z)
        out["Ex_" + tag] = E[0]
        out["Ey_" + tag] = E[1]
        out["Ez_" + tag] = E[2]
        out["Hx_" + tag] = H[0]
        out["Hy_" + tag] = H[1]
        out["Hz_" + tag] = H[2]
        ph = complex(math.cos(2 * math.pi * z / lam), math.sin(2 * math.pi * z / lam))
        out["Ey_phase_err_" + tag] = abs(E[1] - ph)
        out["|E|_" + tag] = abs(E[0]) + abs(E[1]) + abs(E[2])
    # transverse invariance: the field must not depend on x or y here
    E0, _ = S.GetFields(0.0, 0.0, 0.25)
    E1, _ = S.GetFields(0.31, -0.22, 0.25)
    out["xy_invariance"] = max(abs(a - b) for a, b in zip(E0, E1))
    return out


@test
def fields_patterned():
    """Fields inside and around the patterned layer, at fixed sample points."""
    S = rod_sim(120)
    out = {}
    pts = ((0.0, 0.0, 0.075), (0.35, 0.0, 0.075), (0.2, 0.2, 0.075),
           (0.0, 0.0, -0.3), (0.0, 0.0, 0.45))
    for i, (x, y, z) in enumerate(pts):
        E, H = S.GetFields(x, y, z)
        for c, v in zip("xyz", E):
            out["E%s_p%d" % (c, i)] = v
        for c, v in zip("xyz", H):
            out["H%s_p%d" % (c, i)] = v
    return out


@test
def fields_lattice_periodicity():
    """The field must be periodic under a lattice translation at normal incidence."""
    S = rod_sim(120)
    E0, H0 = S.GetFields(0.13, -0.07, 0.075)
    E1, H1 = S.GetFields(1.13, -0.07, 0.075)
    E2, H2 = S.GetFields(0.13, 0.93, 0.075)
    return {"dx_E": max(abs(a - b) for a, b in zip(E0, E1)),
            "dy_E": max(abs(a - b) for a, b in zip(E0, E2)),
            "dx_H": max(abs(a - b) for a, b in zip(H0, H1)),
            "Ey_ref": E0[1], "Hx_ref": H0[0]}


@test
def fields_c4_symmetry():
    """The rod has C4v symmetry; an s-polarised normal wave maps accordingly.

    Rotating the sample point by 90 degrees must rotate the field vector by 90
    degrees, i.e. E_y(x,y) = -E_x(-y,x) for the in-plane components.
    """
    S = rod_sim(120)
    Ea, _ = S.GetFields(0.22, 0.07, 0.075)
    Eb, _ = S.GetFields(-0.07, 0.22, 0.075)
    return {"Ex_a": Ea[0], "Ey_a": Ea[1], "Ex_b": Eb[0], "Ey_b": Eb[1],
            "rot_err_x": abs(Eb[0] + Ea[1]), "rot_err_y": abs(Eb[1] - Ea[0]),
            "Ez_a": Ea[2], "Ez_b": Eb[2]}


@test
def fields_on_grid():
    """GetFieldsOnGrid must agree with GetFields at the same sample points.

    The grid is produced by an inverse FFT of the mode amplitudes, so it only
    reproduces the pointwise field once the grid resolves every harmonic in the
    basis.  NumBasis 80 reaches |m| = 5, so n must be at least 11; below that
    the two disagree for a perfectly good reason and the check means nothing.
    """
    S = rod_sim(80)
    n = 16
    E, H = S.GetFieldsOnGrid(z=0.075, NumSamples=(n, n), Format="Array")
    out = {"nrows": float(len(E)), "ncols": float(len(E[0]))}
    # The row/column convention is not the same in every build, so measure both
    # readings rather than assuming one: whichever is at rounding level is the
    # convention that build uses, and a build where neither is small is broken.
    worst = 0.0
    worst_t = 0.0
    for i in (0, 3, 11):
        for j in (0, 2, 15):
            a = i / float(n)
            b = j / float(n)
            got = E[i][j]
            Eij, _ = S.GetFields(a, b, 0.075)
            Eji, _ = S.GetFields(b, a, 0.075)
            worst = max(worst, max(abs(u - v) for u, v in zip(got, Eij)))
            worst_t = max(worst_t, max(abs(u - v) for u, v in zip(got, Eji)))
    out["grid_vs_point_max_err"] = worst
    out["grid_vs_point_transposed_err"] = worst_t
    out["either_convention_ok"] = 1.0 if min(worst, worst_t) < 1e-9 else 0.0
    out["first_index_is_x"] = 1.0 if worst < worst_t else 0.0
    out["E00_x"] = E[0][0][0]
    out["E00_y"] = E[0][0][1]
    out["E00_z"] = E[0][0][2]
    out["H00_x"] = H[0][0][0]
    return out


@test
def epsilon_readout():
    """GetEpsilon must report the material actually placed at each point."""
    S = rod_sim(400)
    out = {}
    for tag, (x, y) in (("centre", (0.0, 0.0)), ("edge_in", (0.25, 0.0)),
                        ("corner", (0.48, 0.48)), ("mid", (0.0, 0.35))):
        out["eps_" + tag] = S.GetEpsilon(x, y, 0.075)
    out["eps_above"] = S.GetEpsilon(0.0, 0.0, -0.2)
    out["eps_below"] = S.GetEpsilon(0.0, 0.0, 0.4)
    return out


# --------------------------------------------------------------------------
# 8. amplitudes and the S-matrix
# --------------------------------------------------------------------------

@test
def amplitudes_incident():
    """In the incidence layer, the forward order-0 amplitude is the input.

    The excitation sets a unit s-amplitude in the specular order; GetAmplitudes
    on the incidence layer must show exactly that, with every other forward
    amplitude zero.
    """
    S = rod_sim(120)
    fwd, back = S.GetAmplitudes(Layer="top", zOffset=0.0)
    n2 = len(fwd)
    out = {"len_fwd": float(n2), "len_back": float(len(back))}
    out["fwd_max"] = max(abs(a) for a in fwd)
    # order 0 sits first in each of the two polarisation halves
    out["fwd_0"] = fwd[0]
    out["fwd_half"] = fwd[n2 // 2]
    others = [abs(a) for k, a in enumerate(fwd) if k not in (0, n2 // 2)]
    out["fwd_other_max"] = max(others)
    out["back_max"] = max(abs(a) for a in back)
    return out


@test
def amplitudes_transmitted():
    """Transmitted order-0 amplitude must reproduce the transmitted power."""
    S = rod_sim(120)
    inc, _ = S.GetPowerFlux(Layer="top")
    tr, _ = S.GetPowerFlux(Layer="bot")
    fwd, back = S.GetAmplitudes(Layer="bot", zOffset=0.0)
    n2 = len(fwd)
    p = sum(abs(a) ** 2 for a in fwd)
    return {"T_from_flux": tr.real / inc.real,
            "sum|a|^2": p,
            "back_max": max(abs(a) for a in back),
            "a0": fwd[0], "a_half": fwd[n2 // 2]}


@test
def smatrix_determinant():
    """GetSMatrixDeterminant returns (mantissa, base, exponent).

    The determinant of the full 4n x 4n system matrix underflows long before a
    useful basis size, so the value only carries information at a small n; at
    NumBasis 120 both builds report a mantissa near zero and the comparison
    says nothing.  Small bases are used here so the digits are real.
    """
    out = {}
    for nb in (1, 5, 9):
        S = rod_sim(nb)
        m, base, expo = S.GetSMatrixDeterminant()
        out["mant_nb%d" % nb] = m
        out["base_nb%d" % nb] = float(base)
        out["expo_nb%d" % nb] = float(expo)
        out["log_nb%d" % nb] = (math.log(abs(m)) / math.log(base) + expo
                                if abs(m) > 0 else float("nan"))
    return out


# --------------------------------------------------------------------------
# 9. integrals
# --------------------------------------------------------------------------

@test
def volume_integrals():
    """All four GetLayerVolumeIntegral quantities on the patterned layer."""
    S = rod_sim(120)
    out = {}
    for q in ("U", "E", "H", "e"):
        out["int_" + q] = S.GetLayerVolumeIntegral(Layer="pat", Quantity=q)
    out["int_pat_U"] = out["int_U"]
    out["int_gap_e"] = S.GetLayerVolumeIntegral(Layer="top", Quantity="e")
    return out


@test
def z_integral():
    """GetLayerZIntegral at a few in-plane positions."""
    S = rod_sim(120)
    out = {}
    for tag, xy in (("centre", (0.0, 0.0)), ("outside", (0.45, 0.45))):
        e, h = S.GetLayerZIntegral(Layer="pat", xy=xy)
        for c, v in zip("xyz", e):
            out["E%s_%s" % (c, tag)] = v
        for c, v in zip("xyz", h):
            out["H%s_%s" % (c, tag)] = v
    return out


@test
def stress_tensor():
    """GetStressTensorIntegral above, inside and below the patterned layer."""
    S = rod_sim(120)
    out = {}
    for tag, (layer, z) in (("top", ("top", 0.0)), ("pat0", ("pat", 0.0)),
                            ("pat_mid", ("pat", 0.075)), ("bot", ("bot", 0.0))):
        T = S.GetStressTensorIntegral(Layer=layer, zOffset=z)
        for c, v in zip("xyz", T):
            out["T%s_%s" % (c, tag)] = v
    return out


# --------------------------------------------------------------------------
# 10. excitation variants
# --------------------------------------------------------------------------

POL = {}


def set_exterior(S, exc):
    """Call SetExcitationExterior, spelling the polarisation the way this build wants.

    The pre-cb47b74 binding reads the polarisation with PyString_Check, which in
    Python 3 accepts bytes and rejects str; the ported binding takes either.  The
    excitation itself is identical, so the spelling is normalised here rather
    than losing the comparison over it.
    """
    def spell(kind):
        return dict(Excitations=tuple(
            (g, (p.encode() if kind is bytes else p), re, im) for g, p, re, im in exc))

    if POL.get("kind") is None:
        try:
            S.SetExcitationExterior(**spell(str))
            POL["kind"] = str
            return
        except TypeError:
            POL["kind"] = bytes
    S.SetExcitationExterior(**spell(POL["kind"]))


@test
def excitation_exterior():
    """SetExcitationExterior with the specular s order must equal a planewave.

    Both describe the same incident field, so every derived quantity must agree.
    """
    def run(exterior):
        S = rod_sim(120)
        if exterior:
            set_exterior(S, ((1, 'x', 1.0, 0.0),))
        return flux(S)

    a = run(False)
    b = run(True)
    return {"T_planewave": a[0], "T_exterior": b[0], "diff": abs(a[0] - b[0]),
            "R_planewave": a[1], "R_exterior": b[1]}


@test
def excitation_exterior_orders():
    """Exterior excitation in a non-specular order, and in y polarisation."""
    out = {}
    for tag, exc in (("g1x", ((2, 'x', 1.0, 0.0),)),
                     ("g0y", ((1, 'y', 1.0, 0.0),)),
                     ("mix", ((1, 'x', 0.6, 0.0), (1, 'y', 0.0, 0.8)))):
        S = rod_sim(120)
        set_exterior(S, exc)
        f, b = S.GetPowerFlux(Layer="top")
        t, _ = S.GetPowerFlux(Layer="bot")
        out["fwd_" + tag] = f.real
        out["bwd_" + tag] = b.real
        out["trans_" + tag] = t.real
        if abs(f.real) > 1e-30:
            out["ratio_" + tag] = (t.real + abs(b.real)) / f.real
    return out


@test
def complex_frequency():
    """A complex frequency must run and must not silently drop the imaginary part.

    SetFrequency takes one argument parsed as a Python complex, not a pair of
    reals, so the imaginary part travels inside the number.  S4 warns about a
    positive imaginary part, hence the negative sign.
    """
    out = {}
    S = rod_sim(80)
    S.SetFrequency(1 / 1.6)
    out["T_real"] = flux(S)[0]

    for tag, im in (("small", -1e-3), ("big", -1e-2)):
        S = rod_sim(80)
        S.SetFrequency(complex(1 / 1.6, im))
        f, b = S.GetPowerFlux(Layer="top")
        t, _ = S.GetPowerFlux(Layer="bot")
        out["T_" + tag] = t.real / f.real
        out["RpT_" + tag] = t.real / f.real + abs(b.real) / f.real
    out["shift_small"] = abs(out["T_small"] - out["T_real"])
    out["shift_big"] = abs(out["T_big"] - out["T_real"])
    return out


@test
def absorbing_material():
    """A lossy layer: R+T must fall below 1 by the absorbed fraction."""
    S = New(Lattice=((1, 0), (0, 1)), NumBasis=40)
    S.SetMaterial(Name="air", Epsilon=1.0)
    S.SetMaterial(Name="lossy", Epsilon=complex(4.0, 0.25))
    S.AddLayer(Name="top", Thickness=0, Material="air")
    S.AddLayer(Name="mid", Thickness=0.6, Material="lossy")
    S.AddLayer(Name="bot", Thickness=0, Material="air")
    S.SetExcitationPlanewave(IncidenceAngles=(0, 0), sAmplitude=1, pAmplitude=0)
    S.SetFrequency(1 / 1.31)
    T, R = flux(S)
    return {"T": T, "R": R, "RpT": R + T, "absorbed": 1 - (R + T)}


# --------------------------------------------------------------------------
# 11. 1D gratings -- master's two most recent commits are 1D patterning fixes
# --------------------------------------------------------------------------

def _grating_1d(nb, basis="Normal", period=1.0, width=0.5, thick=0.15, lam=1.6,
                theta=0.0):
    S = New(Lattice=period, NumBasis=nb)
    if basis is not None:
        S.SetOptions(PolarizationDecomposition=True, PolarizationBasis=basis)
    S.SetMaterial(Name="air", Epsilon=1.0)
    S.SetMaterial(Name="rod", Epsilon=3.5 ** 2)
    S.AddLayer(Name="top", Thickness=0, Material="air")
    S.AddLayer(Name="pat", Thickness=thick, Material="air")
    S.SetRegionRectangle(Layer="pat", Material="rod", Center=(0, 0), Angle=0,
                         Halfwidths=(width / 2.0, 0))
    S.AddLayer(Name="bot", Thickness=0, Material="air")
    S.SetExcitationPlanewave(IncidenceAngles=(theta, 0), sAmplitude=1, pAmplitude=0)
    S.SetFrequency(1 / lam)
    return S


@test
def grating_1d_basic():
    """A 1D lamellar grating: transmission and energy conservation."""
    out = {}
    for nb in (11, 31, 51, 101):
        S = _grating_1d(nb)
        T, R = flux(S)
        out["T_nb%d" % nb] = T
        out["RpT_nb%d" % nb] = R + T
    return out


@test
def grating_1d_vs_2d():
    """A 1D grating must agree with the same stripe in a 2D lattice.

    The 2D run uses a rectangle spanning the whole unit cell in y, which is the
    same structure; the two solvers differ only in how they enumerate G.
    """
    S1 = _grating_1d(101)
    T1, R1 = flux(S1)

    S2 = New(Lattice=((1, 0), (0, 1)), NumBasis=101)
    S2.SetOptions(PolarizationDecomposition=True, PolarizationBasis="Normal")
    S2.SetMaterial(Name="air", Epsilon=1.0)
    S2.SetMaterial(Name="rod", Epsilon=3.5 ** 2)
    S2.AddLayer(Name="top", Thickness=0, Material="air")
    S2.AddLayer(Name="pat", Thickness=0.15, Material="air")
    S2.SetRegionRectangle(Layer="pat", Material="rod", Center=(0, 0), Angle=0,
                          Halfwidths=(0.25, 0.5))
    S2.AddLayer(Name="bot", Thickness=0, Material="air")
    S2.SetExcitationPlanewave(IncidenceAngles=(0, 0), sAmplitude=1, pAmplitude=0)
    S2.SetFrequency(1 / 1.6)
    T2, R2 = flux(S2)
    return {"T_1d": T1, "T_2d": T2, "diff": abs(T1 - T2),
            "RpT_1d": R1 + T1, "RpT_2d": R2 + T2}


@test
def grating_1d_uniform_limit():
    """A 1D 'grating' of pure air must reproduce free space exactly."""
    S = New(Lattice=1.0, NumBasis=21)
    S.SetMaterial(Name="air", Epsilon=1.0)
    S.AddLayer(Name="top", Thickness=0, Material="air")
    S.AddLayer(Name="mid", Thickness=0.5, Material="air")
    S.AddLayer(Name="bot", Thickness=0, Material="air")
    S.SetExcitationPlanewave(IncidenceAngles=(0, 0), sAmplitude=1, pAmplitude=0)
    S.SetFrequency(1 / 1.31)
    T, R = flux(S)

    S2 = New(Lattice=1.0, NumBasis=21)
    S2.SetMaterial(Name="air", Epsilon=1.0)
    S2.SetMaterial(Name="glass", Epsilon=2.25)
    S2.AddLayer(Name="top", Thickness=0, Material="air")
    S2.AddLayer(Name="mid", Thickness=0.5, Material="glass")
    S2.AddLayer(Name="bot", Thickness=0, Material="air")
    S2.SetExcitationPlanewave(IncidenceAngles=(0, 0), sAmplitude=1, pAmplitude=0)
    S2.SetFrequency(1 / 1.5)
    Tg, Rg = flux(S2)
    return {"T_vacuum": T, "R_vacuum": R, "T_halfwave": Tg,
            "halfwave_err": abs(Tg - 1.0)}


@test
def grating_1d_orders():
    """1D grating with several open orders: per-order power and the sum."""
    S = _grating_1d(101, lam=0.7)
    G = S.GetBasisSet()
    inc, _ = S.GetPowerFlux(Layer="top")
    by = S.GetPowerFluxByOrder(Layer="bot")
    out = {"n_orders": float(len(by)), "n_basis": float(len(G)),
           "inc": inc.real}
    idx = {}
    for i, g in enumerate(G):
        m = int(g[0]) if not hasattr(g, "__len__") else int(g[0])
        idx.setdefault(m, i)
    for m in (-1, 0, 1):
        if m in idx:
            out["T_%d" % m] = by[idx[m]][0].real / inc.real
    out["T_total"] = sum(o[0].real for o in by) / inc.real
    if "T_-1" in out and "T_1" in out:
        out["sym_pm1"] = abs(out["T_-1"] - out["T_1"])
    return out


@test
def grating_1d_oblique():
    """1D grating off normal: energy conservation and the +/- order split."""
    out = {}
    for theta in (0.0, 20.0, 40.0):
        S = _grating_1d(101, theta=theta)
        T, R = flux(S)
        out["T_th%g" % theta] = T
        out["RpT_th%g" % theta] = R + T
    return out


@test
def grating_1d_fields():
    """Fields from a 1D solve: periodicity in x, invariance in y."""
    S = _grating_1d(51)
    E0, H0 = S.GetFields(0.13, 0.0, 0.075)
    E1, _ = S.GetFields(1.13, 0.0, 0.075)
    E2, _ = S.GetFields(0.13, 0.41, 0.075)
    out = {"period_err": max(abs(a - b) for a, b in zip(E0, E1)),
           "y_invariance": max(abs(a - b) for a, b in zip(E0, E2))}
    for c, v in zip("xyz", E0):
        out["E" + c] = v
    for c, v in zip("xyz", H0):
        out["H" + c] = v
    return out


# --------------------------------------------------------------------------
# 12. misc plumbing
# --------------------------------------------------------------------------

@test
def clone_matches():
    """Clone must produce a simulation that answers identically.

    The clone is asked for its numbers first and the original second, so a
    clone that shares state with its source shows up as damage to the source
    rather than only as a wrong answer from the copy.
    """
    S = rod_sim(80)
    T0, R0 = flux(S)
    C = Sim(S.raw.Clone())
    inc, back = C.GetPowerFlux(Layer="top")
    fwd, _ = C.GetPowerFlux(Layer="bot")
    T1 = fwd.real / inc.real
    R1 = abs(back.real) / inc.real
    T2, R2 = flux(S)
    return {"T_orig": T0, "T_clone": T1, "diff": abs(T0 - T1),
            "R_orig": R0, "R_clone": R1,
            "T_orig_after_clone": T2, "source_damage": abs(T0 - T2)}


@test
def clone_full_structure():
    """Clone a simulation that uses every owning field: names, patterns,
    polygon vertices, a tensor material, and an exterior excitation."""
    S = New(Lattice=((1, 0), (0, 1)), NumBasis=60)
    S.SetOptions(PolarizationDecomposition=True, PolarizationBasis="Normal")
    S.SetMaterial(Name="air", Epsilon=1.0)
    S.SetMaterial(Name="lossy", Epsilon=complex(6.0, 0.2))
    S.SetMaterial(Name="aniso", Epsilon=((9.0, complex(0.5, 0.3), 0),
                                         (complex(-0.5, -0.3), 12.25, 0),
                                         (0, 0, 9.0)))
    S.AddLayer(Name="top", Thickness=0, Material="air")
    S.AddLayer(Name="pat", Thickness=0.15, Material="air")
    S.SetRegionPolygon(Layer="pat", Material="aniso", Center=(0.05, -0.03),
                       Angle=20,
                       Vertices=((-0.22, -0.12), (0.24, -0.16),
                                 (0.19, 0.21), (-0.14, 0.23)))
    S.SetRegionCircle(Layer="pat", Material="lossy", Center=(0.05, -0.03),
                      Radius=0.05)
    S.AddLayerCopy(Name="pat2", Thickness=0.15, Layer="pat")
    S.AddLayer(Name="bot", Thickness=0, Material="air")
    S.SetExcitationPlanewave(IncidenceAngles=(25, 40), sAmplitude=0.6,
                             pAmplitude=0.8)
    S.SetFrequency(1 / 1.6)

    T0, R0 = flux(S)
    eps0 = S.GetEpsilon(0.05, -0.03, 0.075)
    C = Sim(S.raw.Clone())
    inc, back = C.GetPowerFlux(Layer="top")
    fwd, _ = C.GetPowerFlux(Layer="bot")
    T1 = fwd.real / inc.real
    epsc = C.GetEpsilon(0.05, -0.03, 0.075)
    T2, _ = flux(S)
    return {"T_orig": T0, "T_clone": T1, "diff": abs(T0 - T1),
            "eps_orig": eps0, "eps_clone": epsc,
            "eps_diff": abs(eps0 - epsc),
            "T_orig_after_clone": T2, "source_damage": abs(T0 - T2)}


@test
def clone_is_independent():
    """Editing the clone must not reach back into the original."""
    S = rod_sim(60)
    T0 = flux(S)[0]
    C = Sim(S.raw.Clone())
    C.SetLayerThickness(Layer="pat", Thickness=0.40)
    C.SetMaterial(Name="rod", Epsilon=2.0)
    inc, _ = C.GetPowerFlux(Layer="top")
    fwd, _ = C.GetPowerFlux(Layer="bot")
    Tc = fwd.real / inc.real
    T1 = flux(S)[0]
    return {"T_orig": T0, "T_orig_after_edit": T1, "source_damage": abs(T0 - T1),
            "T_clone_edited": Tc, "clone_actually_changed": abs(Tc - T0)}


@test
def module_entry_points():
    """The module-level functions, which the sweep had never called.

    Two of the three did not work.  SolveInParallel parsed nothing and returned
    None for any arguments or none at all, so a caller was told the layers had
    been solved when nothing had happened.  NewSpectrumSampler's format string
    wrote the optional-argument separator six times, which CPython rejects, so
    it raised SystemError the moment any optional argument was given and only
    worked in its most degenerate form.

    Solving a layer ahead of time must not change what comes out of it, so the
    check is that the transmission is identical either way.
    """
    out = {}

    def rod():
        S = New(Lattice=((1, 0), (0, 1)), NumBasis=40)
        S.SetMaterial(Name="a", Epsilon=1.0)
        S.SetMaterial(Name="r", Epsilon=6.0)
        S.AddLayer(Name="t", Thickness=0, Material="a")
        S.AddLayer(Name="L", Thickness=0.2, Material="a")
        S.SetRegionCircle(Layer="L", Material="r", Center=(0, 0), Radius=0.25)
        S.AddLayer(Name="b", Thickness=0, Material="a")
        S.SetExcitationPlanewave(IncidenceAngles=(0, 0), sAmplitude=1,
                                 pAmplitude=0)
        S.SetFrequency(1 / 1.6)
        return S

    import s4compat
    a, b = rod(), rod()
    try:
        s4compat.S4.SolveInParallel(
            **{("S4_Layer" if s4compat.PREFIXED else "Layer"): "L",
               "Simulations": [a.raw, b.raw]})
        out["solve_in_parallel_ran"] = 1.0
    except Exception:
        out["solve_in_parallel_ran"] = 0.0
    def T(S):
        return flux(S, layer="b", top="t")[0]
    out["T_presolved"] = T(a)
    out["T_plain"] = T(rod())
    out["presolve_changes_nothing"] = abs(out["T_presolved"] - out["T_plain"])

    # it must object to arguments it cannot use, rather than accepting anything
    def rejects(**kw):
        try:
            s4compat.S4.SolveInParallel(**kw)
            return 0.0
        except Exception:
            return 1.0
    lay = "S4_Layer" if s4compat.PREFIXED else "Layer"
    out["rejects_no_args"] = rejects()
    out["rejects_non_sim"] = rejects(**{lay: "L", "Simulations": [a.raw, "x"]})
    out["rejects_bad_layer"] = rejects(**{lay: "nope", "Simulations": [a.raw]})
    out["rejects_non_sequence"] = rejects(**{lay: "L", "Simulations": 42})

    # SpectrumSampler must accept its optional arguments
    for tag, kw in (("required", dict(FreqStart=0.5, FreqEnd=0.6)),
                    ("one_opt", dict(FreqStart=0.5, FreqEnd=0.6,
                                     InitialNumPoints=7)),
                    ("all_opt", dict(FreqStart=0.5, FreqEnd=0.6,
                                     InitialNumPoints=7, RangeThreshold=2e-3,
                                     MaxBend=5.0, MinimumSpacing=1e-5,
                                     Parallelize=False))):
        try:
            sp = s4compat.S4.NewSpectrumSampler(**kw)
            out["sampler_" + tag] = 1.0
            out["sampler_freq_" + tag] = float(sp.GetFrequency())
        except Exception:
            out["sampler_" + tag] = 0.0
            out["sampler_freq_" + tag] = 0.0

    # the one that already worked
    try:
        it = s4compat.S4.NewInterpolator(
            Type="linear", Table=((0.0, (1.0, 2.0)), (1.0, (3.0, 4.0))))
        g = it.Get(0.5)
        out["interpolator_ok"] = 1.0
        out["interpolator_0"] = float(g[0])
        out["interpolator_1"] = float(g[1])
    except Exception:
        out["interpolator_ok"] = 0.0
    return out


@test
def periodic_image_overlap():
    """A region may leave the cell, but it must not meet its own repeats.

    The Fourier transform is taken over the shape as given, which is right for a
    periodic structure only while the shape and its repeats are disjoint.  Once
    they overlap, the shared area is counted once per copy and the structure
    solved is not the one described -- a circle of radius 0.7 in a 1x1 cell of
    eps 9 on eps 1 reported a mean permittivity of 13.3, when filling the whole
    cell with the rod material cannot exceed 9.

    Touching is not overlapping, and the distinction has to be kept: a rectangle
    spanning the cell exactly is how a 1D grating is written inside a 2D lattice,
    and it meets its neighbour along a line.
    """
    out = {}

    def accepts_circle(rad):
        return _accepts(lambda S: S.SetRegionCircle(
            Layer="pat", Material="e9", Center=(0, 0), Radius=rad))

    def accepts_rect(hx, hy):
        return _accepts(lambda S: S.SetRegionRectangle(
            Layer="pat", Material="e9", Center=(0, 0), Angle=0,
            Halfwidths=(hx, hy)))

    # a circle in a 1x1 cell touches its image at exactly 0.5
    for rad in (0.30, 0.49, 0.4999, 0.50):
        out["circle_%g" % rad] = accepts_circle(rad)          # must be accepted
    for rad in (0.5001, 0.55, 0.70):
        out["circle_over_%g" % rad] = accepts_circle(rad)     # must be refused

    out["rect_fills_y"] = accepts_rect(0.25, 0.50)            # how a 1D grating is written in a 2D lattice
    out["rect_fills_cell"] = accepts_rect(0.50, 0.50)
    out["rect_inside"] = accepts_rect(0.25, 0.30)
    out["rect_over"] = accepts_rect(0.25, 0.51)               # must be refused

    out["touching_kept"] = float(sum(
        1 for k, v in out.items()
        if not k.startswith(("circle_over", "rect_over")) and v == 1.0))
    out["overlaps_caught"] = float(sum(
        1 for k, v in out.items()
        if k.startswith(("circle_over", "rect_over")) and v == 0.0))
    out["n_touching"] = float(sum(
        1 for k in out if k.startswith(("circle_0", "rect_")) and
        not k.startswith("rect_over")))
    out["n_overlaps"] = 4.0
    return out


def _mean_eps_1d(place, n=256, nb=21, thick=0.35):
    """Accepted flag and mean permittivity for a 1D pattern.

    Returns (1.0, mean) if S4 solves it, (0.0, nan) if it refuses.
    """
    S = New(Lattice=1.0, NumBasis=nb)
    S.SetMaterial(Name="bg", Epsilon=1.0)
    S.SetMaterial(Name="e4", Epsilon=4.0)
    S.SetMaterial(Name="e9", Epsilon=9.0)
    S.AddLayer(Name="top", Thickness=0, Material="bg")
    S.AddLayer(Name="pat", Thickness=thick, Material="bg")
    place(S)
    S.AddLayer(Name="bot", Thickness=0, Material="bg")
    S.SetExcitationPlanewave(IncidenceAngles=(0, 0), sAmplitude=1, pAmplitude=0)
    S.SetFrequency(1 / 1.6)
    try:
        S.GetPowerFlux(Layer="top")
    except RuntimeError:
        return 0.0, float("nan")
    z = thick / 2
    tot = sum(S.GetEpsilon((i + 0.5) / n - 0.5, 0.0, z).real for i in range(n))
    return 1.0, tot / n


@test
def reciprocal_lattice():
    """GetReciprocalLattice on square and oblique lattices."""
    out = {}
    S = New(Lattice=((1, 0), (0, 1)), NumBasis=10)
    (a, b), (c, d) = S.GetReciprocalLattice()
    out["sq_00"], out["sq_01"] = a, b
    out["sq_10"], out["sq_11"] = c, d
    S = New(Lattice=((1.0, 0.0), (0.5, 0.8660254037844386)), NumBasis=10)
    (a, b), (c, d) = S.GetReciprocalLattice()
    out["hex_00"], out["hex_01"] = a, b
    out["hex_10"], out["hex_11"] = c, d
    return out


@test
def hex_lattice():
    """A triangular lattice of rods -- non-orthogonal basis vectors."""
    S = New(Lattice=((1.0, 0.0), (0.5, 0.8660254037844386)), NumBasis=91)
    S.SetOptions(PolarizationDecomposition=True, PolarizationBasis="Normal")
    S.SetMaterial(Name="air", Epsilon=1.0)
    S.SetMaterial(Name="rod", Epsilon=3.5 ** 2)
    S.AddLayer(Name="top", Thickness=0, Material="air")
    S.AddLayer(Name="pat", Thickness=0.15, Material="air")
    S.SetRegionCircle(Layer="pat", Material="rod", Center=(0, 0), Radius=0.25)
    S.AddLayer(Name="bot", Thickness=0, Material="air")
    S.SetExcitationPlanewave(IncidenceAngles=(0, 0), sAmplitude=1, pAmplitude=0)
    S.SetFrequency(1 / 1.6)
    T, R = flux(S)
    return {"T": T, "RpT": R + T}


@test
def uniform_slab_fresnel():
    """Uniform slabs against the analytic Fresnel result."""
    out = {}
    for tag, (n, d, lam) in (("halfwave", (1.5, 0.5, 1.5)),
                             ("a", (2.0, 0.37, 1.31)),
                             ("b", (1.5, 0.23, 0.97))):
        S = New(Lattice=((1, 0), (0, 1)), NumBasis=20)
        S.SetMaterial(Name="air", Epsilon=1.0)
        S.SetMaterial(Name="m", Epsilon=n * n)
        S.AddLayer(Name="top", Thickness=0, Material="air")
        S.AddLayer(Name="mid", Thickness=d, Material="m")
        S.AddLayer(Name="bot", Thickness=0, Material="air")
        S.SetExcitationPlanewave(IncidenceAngles=(0, 0), sAmplitude=1, pAmplitude=0)
        S.SetFrequency(1 / lam)
        T, R = flux(S)
        k = 2 * math.pi / lam
        r = (1 - n) / (1 + n)
        e = complex(math.cos(2 * n * k * d), math.sin(2 * n * k * d))
        rtot = r * (1 - e) / (1 - r * r * e)
        out["T_" + tag] = T
        out["Rfresnel_" + tag] = abs(rtot) ** 2
        out["err_" + tag] = abs((1 - T) - abs(rtot) ** 2)
    return out


@test
def output_files():
    """The file-writing entry points must run without error."""
    import os
    import tempfile
    d = tempfile.mkdtemp(prefix="s4probe")
    cwd = os.getcwd()
    os.chdir(d)
    try:
        S = rod_sim(40)
        S.OutputLayerPatternPostscript(Layer="pat", Filename="pat.ps")
        S.OutputLayerPatternRealization(Layer="pat", Nu=8, Nv=8,
                                        Filename="real.txt")
        S.GetPowerFlux(Layer="top")
        S.OutputStructurePOVRay(Filename="struct.pov")
        S.GetFieldsOnGrid(z=0.075, NumSamples=(4, 4), Format="FileWrite",
                          BaseFilename="fld")
        out = {}
        for f in ("pat.ps", "real.txt", "struct.pov", "fld.E", "fld.H"):
            out["size_" + f.replace(".", "_")] = (
                float(os.path.getsize(f)) if os.path.exists(f) else -1.0)
        return out
    finally:
        os.chdir(cwd)


# --------------------------------------------------------------------------
# 13. exact geometry: the DC Fourier coefficient
#
# GetEpsilon returns S4's band-limited reconstruction of the permittivity, so
# its mean over a grid finer than the highest harmonic is the DC coefficient
# exactly -- no truncation error at all.  That makes the mean an exact readout
# of what geometry the solver actually built, independent of the solve.
# --------------------------------------------------------------------------

NB_DC = 25          # |m| <= 3 or so; a 64x64 grid is far past Nyquist
NGRID_DC = 64


def mean_eps(place, nb=NB_DC, n=NGRID_DC, thick=0.35):
    S = New(Lattice=((1, 0), (0, 1)), NumBasis=nb)
    S.SetMaterial(Name="bg", Epsilon=1.0)
    S.SetMaterial(Name="e4", Epsilon=4.0)
    S.SetMaterial(Name="e9", Epsilon=9.0)
    S.AddLayer(Name="top", Thickness=0, Material="bg")
    S.AddLayer(Name="pat", Thickness=thick, Material="bg")
    place(S)
    S.AddLayer(Name="bot", Thickness=0, Material="bg")
    S.SetExcitationPlanewave(IncidenceAngles=(0, 0), sAmplitude=1, pAmplitude=0)
    S.SetFrequency(1 / 1.6)
    z = thick / 2
    try:
        S.GetPowerFlux(Layer="top")
    except RuntimeError:
        return None       # S4 refused this pattern
    tot = 0.0
    for i in range(n):
        x = (i + 0.5) / n - 0.5
        for j in range(n):
            y = (j + 0.5) / n - 0.5
            tot += S.GetEpsilon(x, y, z).real
    return tot / (n * n)


@test
def dc_single_region():
    """One square rod: the mean must be exactly the area-weighted average."""
    hw = 0.17
    got = mean_eps(lambda S: S.SetRegionRectangle(
        Layer="pat", Material="e4", Center=(0, 0), Angle=0, Halfwidths=(hw, hw)))
    want = 1.0 + (4.0 - 1.0) * (2 * hw) ** 2
    return {"mean": got, "expected": want, "err": abs(got - want)}


@test
def dc_disjoint_regions():
    """Two rods that do not touch: both must contribute in full."""
    hw = 0.17

    def place(S):
        S.SetRegionRectangle(Layer="pat", Material="e4", Center=(-0.25, -0.25),
                             Angle=0, Halfwidths=(hw, hw))
        S.SetRegionRectangle(Layer="pat", Material="e9", Center=(0.25, 0.25),
                             Angle=0, Halfwidths=(hw, hw))
    got = mean_eps(place)
    a = (2 * hw) ** 2
    want = 1.0 + 3.0 * a + 8.0 * a
    return {"mean": got, "expected": want, "err": abs(got - want)}


@test
def dc_nested_regions():
    """A ring: the inner region replaces the outer one where they overlap."""
    ho, hi = 0.40, 0.17

    def place(S):
        S.SetRegionRectangle(Layer="pat", Material="e9", Center=(0, 0),
                             Angle=0, Halfwidths=(ho, ho))
        S.SetRegionRectangle(Layer="pat", Material="e4", Center=(0, 0),
                             Angle=0, Halfwidths=(hi, hi))
    got = mean_eps(place)
    ao, ai = (2 * ho) ** 2, (2 * hi) ** 2
    want = 1.0 + 8.0 * ao + (4.0 - 9.0) * ai
    return {"mean": got, "expected": want, "err": abs(got - want)}


@test
def dc_crossing_regions():
    """Two rods whose edges cross -- S4 has no way to express this.

    The containment tree can say "j is inside i" or "j is beside i" and nothing
    else, so a pair that is partly both gets forced into the nesting branch and
    the earlier rod's exposed area silently disappears.  The expected value
    below is the later-region-wins reading, which is what every other region
    pair in S4 obeys.
    """
    hw, off = 0.17, 0.08

    def place(S):
        S.SetRegionRectangle(Layer="pat", Material="e4", Center=(-off, 0),
                             Angle=0, Halfwidths=(hw, hw))
        S.SetRegionRectangle(Layer="pat", Material="e9", Center=(off, 0),
                             Angle=0, Halfwidths=(hw, hw))
    got = mean_eps(place)
    a = (2 * hw) ** 2
    overlap = (2 * hw - 2 * off) * (2 * hw)
    want = 1.0 + 3.0 * (a - overlap) + 8.0 * a
    # what the containment-tree reading produces instead: the eps-4 rod is
    # treated as the parent of the eps-9 rod, so only (9-4) shows through it
    forced_nesting = 1.0 + 3.0 * a + (8.0 - 3.0) * a
    return {"accepted": 0.0 if got is None else 1.0,
            "mean": 0.0 if got is None else got,
            "expected": want,
            "err": 0.0 if got is None else abs(got - want),
            "forced_nesting_value": forced_nesting,
            "err_vs_forced_nesting":
                0.0 if got is None else abs(got - forced_nesting)}


@test
def dc_rotated_offcentre_nesting():
    """Nesting inside a rotated, off-diagonally-centred region.

    shape_contains_point rotates the query point about the region centre.  If
    that rotation reads center[1] where it means center[0], the query is
    displaced by sin(angle)*(cy - cx) -- zero for an unrotated region and zero
    on the diagonal cx == cy, which is why every ordinary test misses it.  Here
    the displacement is 0.212 against a half-width of 0.08, so the inner circle
    is placed outside its parent and the ring becomes an overlap.
    """
    cx, cy, ang = 0.15, -0.15, 45.0
    hx, hy, rad = 0.40, 0.08, 0.04

    def place(S):
        S.SetRegionRectangle(Layer="pat", Material="e9", Center=(cx, cy),
                             Angle=ang, Halfwidths=(hx, hy))
        S.SetRegionCircle(Layer="pat", Material="e4", Center=(cx, cy),
                          Radius=rad)
    got = mean_eps(place)
    ao = 4 * hx * hy
    ai = math.pi * rad * rad
    want = 1.0 + 8.0 * ao + (4.0 - 9.0) * ai
    # what comes out if the circle is not recognised as being inside the
    # rectangle: it is treated as sitting on the background instead
    orphaned = 1.0 + 8.0 * ao + (4.0 - 1.0) * ai
    return {"mean": got, "expected": want, "err": abs(got - want),
            "orphaned_value": orphaned,
            "err_vs_orphaned": abs(got - orphaned)}


@test
def dc_rotated_offcentre_control():
    """The same nesting with the rotation removed -- the control.

    Identical geometry apart from Angle=0, where the suspected displacement is
    exactly zero.  If this passes and the rotated case above does not, the
    rotation is the only difference that can account for it.
    """
    cx, cy = 0.15, -0.15
    hx, hy, rad = 0.40, 0.08, 0.04

    def place(S):
        S.SetRegionRectangle(Layer="pat", Material="e9", Center=(cx, cy),
                             Angle=0, Halfwidths=(hx, hy))
        S.SetRegionCircle(Layer="pat", Material="e4", Center=(cx, cy),
                          Radius=rad)
    got = mean_eps(place)
    ao = 4 * hx * hy
    ai = math.pi * rad * rad
    want = 1.0 + 8.0 * ao + (4.0 - 9.0) * ai
    return {"mean": got, "expected": want, "err": abs(got - want)}


def _accepts(place):
    """1.0 if S4 solves this pattern, 0.0 if it refuses it.

    Every refusal reaches Python as a RuntimeError, so the distinction being
    measured is solved-or-refused, not which message came back.
    """
    S = New(Lattice=((1, 0), (0, 1)), NumBasis=25)
    S.SetMaterial(Name="bg", Epsilon=1.0)
    S.SetMaterial(Name="e4", Epsilon=4.0)
    S.SetMaterial(Name="e9", Epsilon=9.0)
    S.AddLayer(Name="top", Thickness=0, Material="bg")
    S.AddLayer(Name="pat", Thickness=0.35, Material="bg")
    place(S)
    S.AddLayer(Name="bot", Thickness=0, Material="bg")
    S.SetExcitationPlanewave(IncidenceAngles=(0, 0), sAmplitude=1, pAmplitude=0)
    S.SetFrequency(1 / 1.6)
    try:
        S.GetPowerFlux(Layer="top")
        return 1.0
    except RuntimeError:
        return 0.0


def _rect(S, mat, c, ang, hw):
    S.SetRegionRectangle(Layer="pat", Material=mat, Center=c, Angle=ang,
                         Halfwidths=hw)


def _circ(S, mat, c, r):
    S.SetRegionCircle(Layer="pat", Material=mat, Center=c, Radius=r)


def _ell(S, mat, c, ang, hw):
    S.SetRegionEllipse(Layer="pat", Material=mat, Center=c, Angle=ang,
                       Halfwidths=hw)


def _poly(S, mat, c, ang, v):
    S.SetRegionPolygon(Layer="pat", Material=mat, Center=c, Angle=ang,
                       Vertices=v)


_TRI = ((-0.20, -0.14), (0.20, -0.14), (0.0, 0.22))


@test
def crossing_detection():
    """Crossing regions must be refused; nested and disjoint ones accepted.

    Every unordered pair of the four region types appears, once crossing and
    once legal, because the crossing test dispatches on the pair of types and a
    stub in any branch would let that combination through silently.
    """
    cases = {
        # --- crossing: these must be refused -------------------------------
        "cross_rect_rect": lambda S: (
            _rect(S, "e4", (-0.08, 0), 0, (0.17, 0.17)),
            _rect(S, "e9", (0.08, 0), 0, (0.17, 0.17))),
        "cross_rect_circ": lambda S: (
            _rect(S, "e4", (-0.08, 0), 0, (0.17, 0.17)),
            _circ(S, "e9", (0.10, 0), 0.17)),
        "cross_rect_ell": lambda S: (
            _rect(S, "e4", (-0.08, 0), 0, (0.17, 0.17)),
            _ell(S, "e9", (0.10, 0), 20, (0.20, 0.10))),
        "cross_rect_poly": lambda S: (
            _rect(S, "e4", (-0.08, 0), 0, (0.17, 0.17)),
            _poly(S, "e9", (0.10, 0), 0, _TRI)),
        "cross_circ_circ": lambda S: (
            _circ(S, "e4", (-0.10, 0), 0.20),
            _circ(S, "e9", (0.10, 0), 0.20)),
        "cross_circ_ell": lambda S: (
            _circ(S, "e4", (-0.10, 0), 0.20),
            _ell(S, "e9", (0.12, 0), 35, (0.22, 0.09))),
        "cross_circ_poly": lambda S: (
            _circ(S, "e4", (-0.10, 0), 0.20),
            _poly(S, "e9", (0.12, 0), 0, _TRI)),
        "cross_ell_ell": lambda S: (
            _ell(S, "e4", (0, 0), 0, (0.35, 0.10)),
            _ell(S, "e9", (0, 0), 90, (0.35, 0.10))),
        "cross_ell_poly": lambda S: (
            _ell(S, "e4", (-0.10, 0), 0, (0.24, 0.12)),
            _poly(S, "e9", (0.10, 0), 0, _TRI)),
        "cross_poly_poly": lambda S: (
            _poly(S, "e4", (-0.10, 0), 0, _TRI),
            _poly(S, "e9", (0.10, 0), 180, _TRI)),
        "cross_rotated": lambda S: (
            _rect(S, "e4", (0.10, -0.05), 25, (0.30, 0.10)),
            _rect(S, "e9", (0.10, -0.05), 115, (0.30, 0.10))),
        # --- legal: these must keep working --------------------------------
        "ok_single": lambda S: _rect(S, "e4", (0, 0), 0, (0.17, 0.17)),
        "ok_disjoint": lambda S: (
            _rect(S, "e4", (-0.25, -0.25), 0, (0.17, 0.17)),
            _rect(S, "e9", (0.25, 0.25), 0, (0.17, 0.17))),
        "ok_nested_rect": lambda S: (
            _rect(S, "e9", (0, 0), 0, (0.40, 0.40)),
            _rect(S, "e4", (0, 0), 0, (0.17, 0.17))),
        "ok_nested_circ": lambda S: (
            _circ(S, "e9", (0, 0), 0.40),
            _circ(S, "e4", (0, 0), 0.15)),
        "ok_nested_mixed": lambda S: (
            _ell(S, "e9", (0.05, -0.05), 30, (0.40, 0.30)),
            _poly(S, "e4", (0.05, -0.05), 10, _TRI)),
        "ok_nested_rotated": lambda S: (
            _rect(S, "e9", (0.15, -0.15), 45, (0.40, 0.20)),
            _circ(S, "e4", (0.15, -0.15), 0.08)),
        "ok_three_deep": lambda S: (
            _circ(S, "e9", (0, 0), 0.45),
            _circ(S, "e4", (0, 0), 0.30),
            _circ(S, "e9", (0, 0), 0.12)),
        "ok_four_siblings": lambda S: (
            _circ(S, "e4", (-0.25, -0.25), 0.12),
            _circ(S, "e9", (0.25, -0.25), 0.12),
            _circ(S, "e4", (-0.25, 0.25), 0.12),
            _circ(S, "e9", (0.25, 0.25), 0.12)),
        "ok_touching": lambda S: (
            _circ(S, "e4", (-0.20, 0), 0.20),
            _circ(S, "e9", (0.20, 0), 0.20)),
        "ok_nested_touching": lambda S: (
            _circ(S, "e9", (0, 0), 0.40),
            _circ(S, "e4", (0.20, 0), 0.20)),
    }
    return dict((k, _accepts(v)) for k, v in cases.items())


@test
def crossing_detection_thin():
    """How thin an overlap is still caught, and how narrow a gap stays accepted.

    The interior test is where a sliver has to be resolved, and it is searched
    over the intersection of the two bounding boxes -- which for a thin overlap
    IS the sliver, so the grid spacing shrinks with it.  Searching each shape's
    own box instead, as this used to, gives a spacing set by the shape's size:
    two 0.4-wide rectangles sharing a 0.01 strip were accepted.

    The pairs below are built so that the boundaries never cross transversally
    (the edges are collinear, or they meet at a vertex), which is exactly where
    the exact segment test declines to answer and the interior test is the only
    thing standing.

    The last block is the other half of the claim: shapes separated by a gap of
    the same size must keep being accepted.  A detector that refused everything
    would score perfectly on the first block alone.
    """
    out = {}

    def rect(m, c, hw, ang=0):
        return lambda S: S.SetRegionRectangle(Layer="pat", Material=m, Center=c,
                                              Angle=ang, Halfwidths=hw)

    def both(*fs):
        return lambda S: [f(S) for f in fs]

    for d in (1e-2, 1e-3, 1e-4, 1e-5):
        tag = "%g" % d
        # a strip of width d shared along collinear edges
        out["overlap_strip_" + tag] = _accepts(both(
            rect("e4", (-0.20, 0), (0.20, 0.20)),
            rect("e9", (0.20 - d, 0), (0.20, 0.20))))
        # a d x d square shared at one corner
        out["overlap_corner_" + tag] = _accepts(both(
            rect("e4", (-0.20, -0.20), (0.20, 0.20)),
            rect("e9", (0.20 - d, 0.20 - d), (0.20, 0.20))))
        # the same two rectangles pulled apart by d: must stay accepted
        out["gap_strip_" + tag] = _accepts(both(
            rect("e4", (-0.20, 0), (0.20, 0.20)),
            rect("e9", (0.20 + d, 0), (0.20, 0.20))))

    # rotated squares meeting face to face, overlapping by d and short of it by d
    for d in (1e-2, 1e-3):
        h_in = (0.30 - d) / 2 / math.sqrt(2)
        h_out = (0.30 + d) / 2 / math.sqrt(2)
        tag = "%g" % d
        out["overlap_rot45_" + tag] = _accepts(both(
            rect("e4", (-h_in, -h_in), (0.15, 0.15), 45),
            rect("e9", (h_in, h_in), (0.15, 0.15), 45)))
        out["gap_rot45_" + tag] = _accepts(both(
            rect("e4", (-h_out, -h_out), (0.15, 0.15), 45),
            rect("e9", (h_out, h_out), (0.15, 0.15), 45)))

    out["n_overlaps_caught"] = float(
        sum(1 for k, v in out.items() if k.startswith("overlap_") and v == 0.0))
    out["n_overlaps"] = float(
        sum(1 for k in out if k.startswith("overlap_")))
    out["n_gaps_kept"] = float(
        sum(1 for k, v in out.items() if k.startswith("gap_") and v == 1.0))
    out["n_gaps"] = float(sum(1 for k in out if k.startswith("gap_")))
    return out


@test
def crossing_abutting_clipper_output():
    """Polygons that abut along a shared edge must be accepted.

    A design with overlapping shapes has to be resolved by a clipper before S4
    will take it, and clipper output abuts: the two sides of a cut share an
    edge, with each side's vertices computed separately.  The orientation
    determinant that should be zero there comes out as a few ulps of either
    sign, so a segment test without a tolerance reads it as a transversal
    crossing and refuses geometry that does not overlap at all.

    The two quadrilaterals below are lifted verbatim from an independent implementation's
    cross-check suite's resolved 'overlap_tilt' case.  Their interiors are
    disjoint -- a 1200x1200 point-in-polygon sweep finds zero shared area --
    and the determinants in question are 4e-18 and 6e-17 against terms of order
    0.35.  The exact mean permittivity is recorded alongside so that accepting
    them is not confused with accepting them wrongly.
    """
    A = ((0.84022028406711302, 1.3747012476770732),
         (0.43357838453075448, 1.2266957002403713),
         (0.63879047052615567, 0.66288012776882643),
         (1.440220284067113, 0.95457672475124755))
    B = ((0.83967222513263784, 1.3750850026743759),
         (1.6588242694216297, 0.8015085663233299),
         (2.0603277748673618, 1.3749149973256241),
         (1.2411757305783702, 1.9484914336766701))
    P = 2.5
    shift = tuple(tuple((x - P / 2, y - P / 2)) for x, y in ())  # noqa: F841

    def place(S):
        for V, mat in ((A, "e4"), (B, "e9")):
            S.SetRegionPolygon(
                Layer="pat", Material=mat, Center=(0, 0), Angle=0,
                Vertices=tuple((x - P / 2, y - P / 2) for x, y in V))

    S = New(Lattice=((P, 0), (0, P)), NumBasis=49)
    S.SetMaterial(Name="bg", Epsilon=1.0)
    S.SetMaterial(Name="e4", Epsilon=4.0)
    S.SetMaterial(Name="e9", Epsilon=9.0)
    S.AddLayer(Name="top", Thickness=0, Material="bg")
    S.AddLayer(Name="pat", Thickness=0.9, Material="bg")
    place(S)
    S.AddLayer(Name="bot", Thickness=0, Material="bg")
    S.SetExcitationPlanewave(IncidenceAngles=(0, 0), sAmplitude=1, pAmplitude=0)
    S.SetFrequency(1 / 1.0)
    try:
        S.GetPowerFlux(Layer="top")
        accepted = 1.0
    except RuntimeError:
        return {"accepted": 0.0, "mean": 0.0, "expected": 0.0, "err": 0.0}

    # exact: shoelace areas of the two disjoint loops
    def shoelace(V):
        n = len(V)
        return abs(sum(V[i][0] * V[(i + 1) % n][1] - V[(i + 1) % n][0] * V[i][1]
                       for i in range(n))) / 2
    want = 1.0 + ((4.0 - 1.0) * shoelace(A) + (9.0 - 1.0) * shoelace(B)) / (P * P)

    n = 96
    tot = 0.0
    for i in range(n):
        x = (i + 0.5) / n * P - P / 2
        for j in range(n):
            y = (j + 0.5) / n * P - P / 2
            tot += S.GetEpsilon(x, y, 0.45).real
    got = tot / (n * n)
    return {"accepted": accepted, "mean": got, "expected": want,
            "err": abs(got - want)}


def _mean_eps_1d(place, n=256, nb=21, thick=0.35):
    """Accepted flag and mean permittivity for a 1D pattern.

    Returns (1.0, mean) if S4 solves it, (0.0, nan) if it refuses.
    """
    S = New(Lattice=1.0, NumBasis=nb)
    S.SetMaterial(Name="bg", Epsilon=1.0)
    S.SetMaterial(Name="e4", Epsilon=4.0)
    S.SetMaterial(Name="e9", Epsilon=9.0)
    S.AddLayer(Name="top", Thickness=0, Material="bg")
    S.AddLayer(Name="pat", Thickness=thick, Material="bg")
    place(S)
    S.AddLayer(Name="bot", Thickness=0, Material="bg")
    S.SetExcitationPlanewave(IncidenceAngles=(0, 0), sAmplitude=1, pAmplitude=0)
    S.SetFrequency(1 / 1.6)
    try:
        S.GetPowerFlux(Layer="top")
    except RuntimeError:
        return 0.0, float("nan")
    z = thick / 2
    tot = sum(S.GetEpsilon((i + 0.5) / n - 0.5, 0.0, z).real for i in range(n))
    return 1.0, tot / n


@test
def crossing_detection_1d():
    """1D regions are intervals, and S4 can only add them to the background.

    There is no nesting in 1D: a region's area is zero, so the containment tree
    leaves every 1D region parented to the background and the transform adds
    each interval's full contribution.  Two intervals that share any length
    therefore double-count it, and that is true of a nested pair as much as a
    partly overlapping one -- so the 2D distinction between "nested, fine" and
    "crossing, not representable" does not exist here.  Both must be refused.

    The exact means below are what the geometry says; the values S4 produced
    when it accepted these patterns are recorded next to them.
    """
    out = {}
    cases = {
        "single": (lambda S: _rect(S, "e4", (0, 0), 0, (0.25, 0)),
                   1.0 + 3.0 * 0.5),
        "disjoint": (lambda S: (_rect(S, "e4", (-0.25, 0), 0, (0.10, 0)),
                                _rect(S, "e9", (0.25, 0), 0, (0.10, 0))),
                     1.0 + 3.0 * 0.2 + 8.0 * 0.2),
        "abutting": (lambda S: (_rect(S, "e9", (-0.20, 0), 0, (0.10, 0)),
                                _rect(S, "e4", (0.0, 0), 0, (0.10, 0))),
                     1.0 + 8.0 * 0.2 + 3.0 * 0.2),
        # a region written inside another: exact answer 5.9, S4 answered 8.3
        "nested": (lambda S: (_rect(S, "e9", (0, 0), 0, (0.40, 0)),
                              _rect(S, "e4", (0, 0), 0, (0.15, 0))),
                   1.0 + 8.0 * 0.8 + (4.0 - 9.0) * 0.3),
        # partly overlapping: exact answer 4.8, S4 answered 5.4
        "overlap": (lambda S: (_rect(S, "e4", (-0.10, 0), 0, (0.20, 0)),
                               _rect(S, "e9", (0.10, 0), 0, (0.20, 0))),
                    1.0 + 3.0 * (0.4 - 0.2) + 8.0 * 0.4),
    }
    for name, (place, want) in cases.items():
        ok, got = _mean_eps_1d(place)
        out["accepted_" + name] = ok
        out["expected_" + name] = want
        out["mean_" + name] = got if ok else 0.0
        out["err_" + name] = abs(got - want) if ok else 0.0
    return out


def _fft_next_fast_size(n):
    """kiss_fft's rule: the next size that factors into 2s, 3s and 5s."""
    while True:
        m = n
        for p in (2, 3, 5):
            while m % p == 0:
                m //= p
        if m <= 1:
            return n
        n += 1


def _discretized(nb, res):
    S = New(Lattice=((1, 0), (0, 1)), NumBasis=nb)
    S.SetOptions(DiscretizedEpsilon=True, DiscretizationResolution=res,
                 PolarizationDecomposition=True, PolarizationBasis="Normal")
    S.SetMaterial(Name="air", Epsilon=1.0)
    S.SetMaterial(Name="rod", Epsilon=complex(3.5 ** 2, 0))
    S.AddLayer(Name="top", Thickness=0, Material="air")
    S.AddLayer(Name="pat", Thickness=0.15, Material="air")
    S.SetRegionCircle(Layer="pat", Material="rod", Center=(0, 0), Radius=0.30)
    S.AddLayer(Name="bot", Thickness=0, Material="air")
    S.SetExcitationPlanewave(IncidenceAngles=(0, 0), sAmplitude=1, pAmplitude=0)
    S.SetFrequency(1 / 1.6)
    return flux(S)


@test
def discretization_grid_parity():
    """Discretised epsilon is wrong whenever the FFT grid comes out odd.

    DiscretizationResolution is a multiplier: the grid is
    fft_next_fast_size(resolution * Gmax).  At NumBasis 400 (Gmax 11) the
    resolutions that land on an odd grid are 4, 11 and 12, and those three are
    exactly the ones whose transmission leaves the sequence the even grids
    trace.  Resolution 4 is not special -- it is just the first multiplier that
    happens to hit an odd size for this Gmax.
    """
    out = {}
    for res in (2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 13):
        T, R = _discretized(400, res)
        ngrid = _fft_next_fast_size(11 * res)
        out["T_res%d" % res] = T
        out["RpT_res%d" % res] = R + T
        out["oldrule_ngrid_res%d" % res] = float(ngrid)
        # the grid the pre-fix rule produced; the parity of this number is
        # what used to decide whether the answer was an outlier
        out["oldrule_odd_res%d" % res] = float(ngrid % 2)
    even = [v for k, v in out.items()
            if k.startswith("T_res") and out["oldrule_odd_res" + k[5:]] == 0.0]
    out["oldrule_even_spread"] = max(even) - min(even)
    odd = [v for k, v in out.items()
           if k.startswith("T_res") and out["oldrule_odd_res" + k[5:]] == 1.0]
    out["oldrule_odd_spread"] = max(odd) - min(odd)
    out["oldrule_odd_vs_even_gap"] = max(odd) - max(even)
    return out


@test
def discretization_grid_overread():
    """Energy conservation breaks exactly where the grid is too small to index.

    The coefficient read is Fto[f[1] + f[0]*ngrid[1]] with f up to 2*Gmax and
    only the negative half wrapped, so it stays in bounds only while
    ngrid > 2*Gmax.  At resolution 2 the grid is fft_next_fast_size(2*Gmax),
    which equals 2*Gmax whenever that already factors into 2s, 3s and 5s -- and
    2 is the documented minimum resolution.  The two configurations below are
    the ones where that happens; they are also the only ones in this file whose
    R+T is not 1 to machine precision.
    """
    out = {}
    for nb, gmax in ((120, 6), (200, 8), (400, 11)):
        for res in (2, 8):
            T, R = _discretized(nb, res)
            ngrid = _fft_next_fast_size(gmax * res)
            tag = "nb%d_res%d" % (nb, res)
            out["T_" + tag] = T
            out["RpT_" + tag] = R + T
            out["energy_err_" + tag] = abs(R + T - 1.0)
            out["oldrule_ngrid_" + tag] = float(ngrid)
            # 1.0 where the pre-fix grid rule let the read index past the
            # end of the FFT buffer; a static prediction, the same for any
            # build, so what moves between builds is the energy error
            out["predicted_overread_" + tag] = 1.0 if ngrid <= 2 * gmax else 0.0
    return out


@test
def normals_rotated_rectangle():
    """The normal field of a rotated rectangle, read through the solver.

    shape_get_normal's RECTANGLE branch writes -.1 where -1. is meant for the
    negative side.  The function normalises its result before returning, so a
    length-0.1 vector comes back unit-length and pointing the same way; this
    test exists to pin the observable consequence, whatever it turns out to be.
    The Normal polarisation basis is the only route by which those normals
    reach a number Python can see.
    """
    out = {}
    for ang in (0.0, 30.0, 90.0):
        S = New(Lattice=((1, 0), (0, 1)), NumBasis=80)
        S.SetOptions(PolarizationDecomposition=True, PolarizationBasis="Normal")
        S.SetMaterial(Name="air", Epsilon=1.0)
        S.SetMaterial(Name="rod", Epsilon=3.5 ** 2)
        S.AddLayer(Name="top", Thickness=0, Material="air")
        S.AddLayer(Name="pat", Thickness=0.15, Material="air")
        S.SetRegionRectangle(Layer="pat", Material="rod", Center=(0.0, 0.0),
                             Angle=ang, Halfwidths=(0.30, 0.12))
        S.AddLayer(Name="bot", Thickness=0, Material="air")
        S.SetExcitationPlanewave(IncidenceAngles=(20.0, 35.0),
                                 sAmplitude=0.6, pAmplitude=0.8)
        S.SetFrequency(1 / 1.6)
        T, R = flux(S)
        out["T_ang%g" % ang] = T
        out["RpT_ang%g" % ang] = R + T
    return out


@test
def two_layer_interface():
    """A bare interface: two layers and nothing between them.

    Every other probe here builds three layers or more, which is why the whole
    sweep missed that SolveAll drops the incident amplitude when there are
    exactly two.  A plane interface is the simplest structure S4 can express and
    the one case with a closed-form answer, so it is worth its own probe.
    """
    import math
    out = {}
    for n2 in (1.5, 3.5):
        for theta in (0, 45):
            for pol in ("s", "p"):
                S = s4compat.New(Lattice=((1, 0), (0, 1)), NumBasis=9)
                S.SetMaterial(Name="a", Epsilon=1.0)
                S.SetMaterial(Name="b", Epsilon=complex(n2) ** 2)
                S.AddLayer(Name="top", Thickness=0, Material="a")
                S.AddLayer(Name="bot", Thickness=0, Material="b")
                s_amp, p_amp = (1, 0) if pol == "s" else (0, 1)
                S.SetExcitationPlanewave(IncidenceAngles=(theta, 0),
                                         sAmplitude=s_amp, pAmplitude=p_amp)
                S.SetFrequency(1 / 1.6)
                inc, back = S.GetPowerFlux(Layer="top")
                fwd, _ = S.GetPowerFlux(Layer="bot")
                i = complex(inc).real
                key = "n%g_th%d_%s" % (n2, theta, pol)
                out["T_" + key] = complex(fwd).real / i
                out["RpT_" + key] = complex(fwd).real / i + abs(complex(back).real) / i
    return out


# --------------------------------------------------------------------------

def encode(v):
    if isinstance(v, complex):
        return {"__c__": [v.real, v.imag]}
    if isinstance(v, (int, float)):
        return float(v)
    raise TypeError("cannot encode %r" % (v,))


def main(argv):
    if len(argv) >= 2 and argv[1] == "list":
        for k in TESTS:
            print(k)
        return 0
    if len(argv) >= 3 and argv[1] == "run":
        name = argv[2]
        try:
            res = TESTS[name]()
            payload = {"status": "ok", "build": s4compat.where(),
                       "style": s4compat.STYLE,
                       "values": dict((k, encode(v)) for k, v in res.items())}
        except Exception:
            payload = {"status": "error", "build": s4compat.where(),
                       "style": s4compat.STYLE,
                       "traceback": traceback.format_exc()}
        sys.stdout.write("@@JSON@@" + json.dumps(payload) + "\n")
        return 0
    sys.stderr.write(__doc__)
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv))
