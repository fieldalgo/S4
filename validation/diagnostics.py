"""Check that S4 refuses geometry it cannot represent, and says why.

Every other script here compares numbers.  This one compares behaviour, because
the defects it covers do not show up as a wrong number in a converged series --
they show up as a plausible number for a structure the user did not describe:

  * two regions whose boundaries cross are not representable by the containment
    tree, and were silently resolved to whichever reading the sort happened to
    pick (2.0880 became 1.9248);
  * a region overlapping a periodic image is counted once per copy, and a rod
    of radius 0.7 in a 1x1 cell read a mean permittivity of 13.3 when the rod
    material itself is 9.

A refusal is only useful if the caller can act on it, so each case also asserts
what the message says.  The region indices are the order the regions were added
to the layer, which is the only numbering the caller ever sees.

The legitimate cases are asserted too: nesting in 2D, and a region that exactly
spans the cell, must keep working or the checks have been made too strict.  1D
is the exception -- it has no containment tree, and nesting there is refused
rather than silently added up.

    python diagnostics.py                    # against whatever S4 imports
    PYTHONPATH=<old build> python diagnostics.py    # for contrast

Exits non-zero if any case fails.  Run against the 2016 oracle and every
refusal case fails -- that is the record of what was fixed, not a bug here.
"""

import sys

import s4compat

LAYER = "patterned"


def base(lattice=((1, 0), (0, 1)), nb=25):
    S = s4compat.New(Lattice=lattice, NumBasis=nb)
    for name, eps in (("bg", 1.0), ("e4", 4.0), ("e9", 9.0)):
        S.SetMaterial(Name=name, Epsilon=eps)
    S.AddLayer(Name="top", Thickness=0, Material="bg")
    S.AddLayer(Name=LAYER, Thickness=0.35, Material="bg")
    return S


def rect(S, mat, center, halfwidths):
    S.SetRegionRectangle(Center=center, Angle=0, Halfwidths=halfwidths,
                         Layer=LAYER, Material=mat)


def circle(S, mat, center, radius):
    S.SetRegionCircle(Center=center, Radius=radius,
                      Layer=LAYER, Material=mat)


def solve(S):
    """Returns None if the structure was accepted, else the message."""
    S.AddLayer(Name="bot", Thickness=0, Material="bg")
    S.SetExcitationPlanewave(IncidenceAngles=(0, 0), sAmplitude=1, pAmplitude=0)
    S.SetFrequency(1 / 1.6)
    try:
        S.GetPowerFlux(Layer="top")
    except Exception as e:                                   # noqa: BLE001
        return str(e)
    return None


# --------------------------------------------------------------------------
# the cases
# --------------------------------------------------------------------------

def case_cross():
    """Region 3 crosses region 2.  Region 1 is disjoint and must not be named."""
    S = base()
    rect(S, "e4", (-0.35, 0), (0.10, 0.10))
    rect(S, "e9", (0.05, 0), (0.17, 0.17))
    rect(S, "e4", (0.20, 0), (0.17, 0.17))
    return solve(S), ["Region 3", LAYER, "crosses"]


def case_self_image():
    """Radius 0.7 in a 1x1 cell: the rod runs into its own repeats."""
    S = base()
    circle(S, "e9", (0, 0), 0.7)
    return solve(S), ["Region 1", LAYER, "its own periodic image"]


def case_other_image():
    """A small region near the edge landing inside a large region's copy.

    Same double-counting, different pair -- the message has to say so rather
    than claim the region met itself.
    """
    S = base()
    circle(S, "e9", (-0.3, 0), 0.05)
    circle(S, "e9", (0, 0), 0.7)
    return solve(S), ["Region 1", "periodic image of region 2"]


def case_nested_ok():
    """Nesting is exactly what the containment tree is for."""
    S = base()
    circle(S, "e9", (0, 0), 0.2)
    rect(S, "e4", (0, 0), (0.05, 0.05))
    return solve(S), None


def case_spanning_ok():
    """A region spanning the cell meets its neighbour along a line only.

    Touching shares no area, so the transform is still right and this must be
    accepted; refusing it would break every grating written the obvious way.
    """
    S = base()
    rect(S, "e9", (0, 0), (0.5, 0.2))
    return solve(S), None


def case_1d_overlap():
    """1D has no containment tree, so it needs its own interval check."""
    S = base(lattice=1.0, nb=9)
    rect(S, "e4", (-0.1, 0), (0.2, 0))
    rect(S, "e9", (0.1, 0), (0.2, 0))
    return solve(S), ["Regions 1 and 2", LAYER, "overlap"]


def case_1d_nested():
    """Nesting is legal in 2D and not in 1D, so this must be refused.

    The 1D transform adds each interval against the layer background and never
    subtracts the enclosing one.  A 0.3 half-width of eps 9 holding a 0.1
    half-width of eps 4 reads a mean permittivity of 6.4 -- 1 + 8*0.6 + 3*0.2,
    the naive superposition -- where the structure described averages 4.8.
    """
    S = base(lattice=1.0, nb=9)
    rect(S, "e9", (0, 0), (0.3, 0))
    rect(S, "e4", (0, 0), (0.1, 0))
    return solve(S), ["Regions 1 and 2", LAYER, "1D patterning"]


def case_1d_across_boundary():
    """Two intervals far apart as written that meet across the cell boundary.

    Comparing centres without reducing modulo the period misses this, and it is
    the natural way to write a grating whose feature straddles the edge.
    """
    S = base(lattice=1.0, nb=9)
    rect(S, "e9", (-0.45, 0), (0.1, 0))
    rect(S, "e4", (0.45, 0), (0.1, 0))
    return solve(S), ["Region 1", "periodic image of region 2"]


def case_1d_wider_than_cell():
    S = base(lattice=1.0, nb=9)
    rect(S, "e9", (0, 0), (0.6, 0))
    return solve(S), ["Region 1", "wider than the lattice period"]


def case_1d_disjoint_ok():
    S = base(lattice=1.0, nb=9)
    rect(S, "e9", (-0.25, 0), (0.15, 0))
    rect(S, "e4", (0.25, 0), (0.15, 0))
    return solve(S), None


def case_1d_full_cell_ok():
    """Exactly as wide as the cell: the copies touch and share nothing."""
    S = base(lattice=1.0, nb=9)
    rect(S, "e9", (0, 0), (0.5, 0))
    return solve(S), None


def case_unknown_material():
    S = base()
    try:
        rect(S, "no_such_material", (0, 0), (0.1, 0.1))
    except Exception as e:                                   # noqa: BLE001
        return str(e), ["no_such_material", "not found"]
    return None, ["no_such_material", "not found"]


def case_clone():
    """A clone must solve to the same number as the simulation it came from.

    Not a message check, and not the crash reproduction either -- a segfault
    would take the rest of this file down with it, so the SIGSEGV case lives in
    probe.py where each probe runs in its own subprocess.  This is the cheap
    always-on guard: a shallow-copied shape array shows up as disagreement long
    before it shows up as a crash.
    """
    S = base()
    circle(S, "e9", (0, 0), 0.2)
    S.AddLayer(Name="bot", Thickness=0, Material="bg")
    S.SetExcitationPlanewave(IncidenceAngles=(0, 0), sAmplitude=1, pAmplitude=0)
    S.SetFrequency(1 / 1.6)
    try:
        T = s4compat.Sim(S.Clone())
        a = complex(S.GetPowerFlux(Layer="bot")[0]).real
        b = complex(T.GetPowerFlux(Layer="bot")[0]).real
    except Exception as e:                                   # noqa: BLE001
        return "raised %s: %s" % (type(e).__name__, e), None
    if abs(a - b) > 1e-12 * max(1.0, abs(a)):
        return "clone disagrees: %.15g vs %.15g" % (a, b), None
    return None, None


def case_both_keyword_spellings():
    """A script written for either release must run against this build.

    Upstream cb47b74 renamed Layer= and Material= to S4_Layer= and
    S4_Material=, so every script written before it raises TypeError after it
    and vice versa, with no way to tell the builds apart from Python except by
    making a call and watching it fail.  Both spellings are accepted here, and
    they have to produce the same number rather than merely not raising.

    This case is deliberately written against the raw module: s4compat exists
    to paper over exactly this, and would hide what is being tested.
    """
    import S4 as raw

    def run(prefix):
        k = lambda name: prefix + name                       # noqa: E731
        S = raw.New(Lattice=((1, 0), (0, 1)), NumBasis=25)
        S.SetMaterial(Name="bg", Epsilon=1.0)
        S.SetMaterial(Name="rod", Epsilon=9.0)
        S.AddLayer(Name="top", Thickness=0, **{k("Material"): "bg"})
        S.AddLayer(Name="pat", Thickness=0.35, **{k("Material"): "bg"})
        S.SetRegionCircle(Center=(0, 0), Radius=0.2,
                          **{k("Layer"): "pat", k("Material"): "rod"})
        S.AddLayer(Name="bot", Thickness=0, **{k("Material"): "bg"})
        S.SetExcitationPlanewave(IncidenceAngles=(0, 0), sAmplitude=1, pAmplitude=0)
        S.SetFrequency(1 / 1.6)
        return complex(S.GetPowerFlux(**{k("Layer"): "bot"})[0]).real

    try:
        bare, prefixed = run(""), run("S4_")
    except Exception as e:                                   # noqa: BLE001
        return "%s: %s" % (type(e).__name__, e), None
    if bare != prefixed:
        return "the two spellings disagree: %.15g vs %.15g" % (bare, prefixed), None
    return None, None


def case_long_name_message():
    """A long, non-ASCII layer name must not eat the diagnostic.

    The message is built in a fixed buffer, so an unbounded name pushes the
    sentence out of it -- and the part that says what is wrong is at the end,
    which leaves the reader holding the name and no complaint.  The cut is also
    byte-wise, so it lands inside a multi-byte character and Python renders the
    remains as U+FFFD.

    The name is an ordinary English one with a tail of code points chosen purely
    for their encoded lengths -- 2, 3 and 4 bytes -- so the boundary logic is
    exercised at every sequence length rather than only the one the author
    happened to type.  The tail is an encoding fixture, not text.
    """
    tail = "".join(chr(c) for c in (0x00e9, 0x20ac, 0x1d11e))    # 2, 3, 4 bytes
    name = "top_patterned_grating_layer_" + tail * 40
    # Built here rather than through base(), which names the layer for the
    # other cases; the name is the whole point of this one.
    S = s4compat.New(Lattice=((1, 0), (0, 1)), NumBasis=25)
    for n, eps in (("bg", 1.0), ("e9", 9.0)):
        S.SetMaterial(Name=n, Epsilon=eps)
    S.AddLayer(Name="top", Thickness=0, Material="bg")
    S.AddLayer(Name=name, Thickness=0.35, Material="bg")
    S.SetRegionCircle(Layer=name, Material="e9", Center=(-0.3, 0), Radius=0.05)
    S.SetRegionCircle(Layer=name, Material="e9", Center=(0, 0), Radius=0.7)
    S.AddLayer(Name="bot", Thickness=0, Material="bg")
    S.SetExcitationPlanewave(IncidenceAngles=(0, 0), sAmplitude=1, pAmplitude=0)
    S.SetFrequency(1 / 1.6)
    try:
        S.GetPowerFlux(Layer="top")
    except Exception as e:                                   # noqa: BLE001
        msg = str(e)
    else:
        return None, ["periodic image of region 2"]
    if "\ufffd" in msg:
        return msg.replace("\ufffd", "<U+FFFD>"), ["no replacement character"]
    return msg, ["Region 1", "periodic image of region 2",
                 "order they were added to the layer"]


def case_set_verbosity():
    """SetVerbosity must set the verbosity, and reject a level it cannot.

    The parse check was written without its negation, so the success path
    returned NULL with no exception set -- CPython turns that into SystemError
    -- and the failure path fell through to use an uninitialised level.  The
    option was therefore unreachable from Python in every build.

    Asserted through the raw module: this is about the binding, and the shim
    would only add a layer between the call and what it is testing.
    """
    import S4 as raw
    S = raw.New(Lattice=((1, 0), (0, 1)), NumBasis=9)
    try:
        if S.SetVerbosity(Level=3) is not None:
            return "SetVerbosity returned something other than None", None
        S.SetVerbosity(0)                                    # positional too
    except Exception as e:                                   # noqa: BLE001
        return "%s: %s" % (type(e).__name__, e), None
    for level in (-1, 10):
        try:
            S.SetVerbosity(Level=level)
        except TypeError:
            continue
        except Exception as e:                               # noqa: BLE001
            return "level %d raised %s" % (level, type(e).__name__), None
        return "level %d was accepted" % level, None
    return None, None


def case_introspection():
    """The build must be able to say what it is and what it was told.

    Upstream exposes no version and no way to read the options back, so a
    caller cannot confirm that SetOptions was understood -- which is how
    LatticeTruncation being silently ignored survived as long as it did.

    __version__ is this fork's own release plus the commit the binary was built
    from, with a .dirty suffix when the working tree had uncommitted changes
    under S4/.  A result that cannot be traced to a build is not reproducible.
    Upstream's release is kept in __upstream__ rather than reported as ours.

    The round trip is the real assertion: what GetOptions reports has to be
    accepted by SetOptions and produce the same reading.  NumBasis against
    NumBasisRequested is the second: circular truncation is the default and it
    returns fewer orders than were asked for, which is invisible otherwise.
    """
    import S4 as raw
    for attr in ("__version__", "__release__", "__build__", "__upstream__"):
        if not getattr(raw, attr, None):
            return "the module exposes no %s" % attr, None
    if raw.__version__ != "%s+%s" % (raw.__release__, raw.__build__):
        return ("__version__ %r is not __release__ plus __build__"
                % raw.__version__), None
    if raw.__build__ == "unknown":
        return "built outside a checkout: __build__ is 'unknown'", None

    def sim(nb=289):
        S = raw.New(Lattice=((1, 0), (0, 1)), NumBasis=nb)
        S.SetMaterial(Name="bg", Epsilon=1.0)
        S.AddLayer(Name="a", Thickness=0, Material="bg")
        S.AddLayer(Name="b", Thickness=0, Material="bg")
        return S

    strip = lambda d: {k: v for k, v in d.items()                # noqa: E731
                       if k not in ("NumBasis", "NumBasisRequested")}
    S = sim(25)
    S.SetOptions(PolarizationDecomposition=True, PolarizationBasis="Normal",
                 DiscretizationResolution=16, Verbosity=2)
    want = strip(S.GetOptions())
    T = sim(25)
    T.SetOptions(**want)
    if strip(T.GetOptions()) != want:
        return "GetOptions does not round-trip through SetOptions", None

    for trunc, exact in (("Circular", False), ("Parallelogramic", True)):
        S = sim()
        S.SetOptions(LatticeTruncation=trunc)
        g = S.GetOptions()
        if g["LatticeTruncation"] != trunc:
            return "asked for %s, got %s" % (trunc, g["LatticeTruncation"]), None
        if g["NumBasis"] != len(S.GetBasisSet()):
            return "NumBasis disagrees with GetBasisSet under %s" % trunc, None
        if exact and g["NumBasis"] != g["NumBasisRequested"]:
            return "%s should keep all 289 orders, kept %d" % (trunc, g["NumBasis"]), None
        if not exact and g["NumBasis"] >= g["NumBasisRequested"]:
            return "circular truncation should drop orders, kept them all", None
    return None, None


def case_docstrings():
    """help() must not name functions that do not exist.

    Every SetRegion* entry carried its pre-rename Lua name, so the one piece of
    documentation available from inside Python pointed at an API that has not
    existed for years.
    """
    import S4 as raw
    S = raw.New(Lattice=((1, 0), (0, 1)), NumBasis=9)
    wrong = []
    for name in sorted(n for n in dir(S) if not n.startswith("_")):
        doc = (getattr(S, name).__doc__ or "").strip()
        if not doc:
            wrong.append("%s has none" % name)
        elif doc.split("(")[0] != name:
            wrong.append("%s says %r" % (name, doc.split("(")[0]))
    if wrong:
        return "docstrings: " + "; ".join(wrong[:4]), None
    return None, None


def case_two_layer_interface():
    """A plane interface must reproduce Fresnel.

    SolveAll folds both boundary columns into the same block row when there are
    exactly two layers -- (nlayers-1)*n4-n2 is n2 for nlayers == 2 -- so the
    second write overwrote the first and the incident amplitude was discarded.
    Air against glass came back with no reflection, no transmission and
    R+T = 0, silently.

    This is asserted against the closed form rather than against another build,
    because both builds have to agree with physics and neither is an authority.
    Both solve paths are checked: ConserveMemory takes SolveInterior, which was
    always right, and its agreement is what localised the defect.
    """
    import cmath
    import S4 as raw

    def s4_interface(n2, theta, pol, conserve):
        S = raw.New(Lattice=((1, 0), (0, 1)), NumBasis=9)
        S.SetOptions(ConserveMemory=conserve)
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
        return complex(fwd).real / i, abs(complex(back).real) / i

    def fresnel(n2, theta, pol):
        th = cmath.pi * theta / 180
        c1, s1 = cmath.cos(th), cmath.sin(th)
        s2 = s1 / n2
        c2 = cmath.sqrt(1 - s2 * s2)
        r = ((c1 - n2 * c2) / (c1 + n2 * c2) if pol == "s"
             else (n2 * c1 - c2) / (n2 * c1 + c2))
        return 1 - abs(r) ** 2, abs(r) ** 2

    worst = 0.0
    for n2 in (1.5, 3.5):
        for theta in (0, 20, 45, 70):
            for pol in ("s", "p"):
                want_T, want_R = fresnel(n2, theta, pol)
                for conserve in (False, True):
                    T, R = s4_interface(n2, theta, pol, conserve)
                    worst = max(worst, abs(T - want_T), abs(R - want_R))
                    if abs(T + R - 1) > 1e-10:
                        return ("n=%g theta=%d %s ConserveMemory=%s gives R+T=%.12g"
                                % (n2, theta, pol, conserve, T + R)), None
    if worst > 1e-12:
        return "worst disagreement with Fresnel is %.3e" % worst, None
    return None, None


def case_duplicate_layer_name():
    """Two layers of the same name both change the answer; only one is reachable.

    Every lookup resolves to the first match, so the second cannot be
    addressed, patterned or measured, while it still sits in the stack and
    moves the result.  Upstream reserved error code 12 for this and left the
    check commented out -- and inside the non-copy branch, where a pair of
    copies would have slipped past it.

    Redefining a layer is what SetLayer is for, so a repeated name from
    AddLayer is a mistake rather than an idiom; case_redefine_layer_ok below
    holds that door open.
    """
    S = base()
    S.AddLayer(Name=LAYER, Thickness=0.5, Material="e9")     # LAYER again
    return solve(S), ["are both named", LAYER]


def case_redefine_layer_ok():
    """SetLayer redefines rather than appends, and must keep working."""
    S = base()
    S.SetLayer(Name=LAYER, Thickness=0.9, Material="e9")
    return solve(S), None


def case_degenerate_lattice():
    """Parallel basis vectors span no cell, so there is no reciprocal lattice.

    S4_Lattice_Reciprocate detects it and returns non-zero; S4_Simulation_New
    ignores the return, and the code that would have reported it is commented
    out there.  What the caller saw instead was LAPACK refusing the layer
    eigensystem on stderr -- not as an exception -- and NaN for every result.
    """
    S = s4compat.New(Lattice=((1, 0), (2, 0)), NumBasis=9)
    S.SetMaterial(Name="bg", Epsilon=1.0)
    S.AddLayer(Name="top", Thickness=0, Material="bg")
    S.AddLayer(Name="bot", Thickness=0, Material="bg")
    S.SetExcitationPlanewave(IncidenceAngles=(0, 0), sAmplitude=1, pAmplitude=0)
    S.SetFrequency(1 / 1.6)
    try:
        S.GetPowerFlux(Layer="top")
    except Exception as e:                                   # noqa: BLE001
        return str(e), ["parallel"]
    return None, ["parallel"]


def case_zero_period_1d():
    """A 1D lattice of period zero is the same degeneracy, spelled shorter."""
    S = s4compat.New(Lattice=0.0, NumBasis=9)
    S.SetMaterial(Name="bg", Epsilon=1.0)
    S.AddLayer(Name="top", Thickness=0, Material="bg")
    S.AddLayer(Name="bot", Thickness=0, Material="bg")
    S.SetExcitationPlanewave(IncidenceAngles=(0, 0), sAmplitude=1, pAmplitude=0)
    S.SetFrequency(1 / 1.6)
    try:
        S.GetPowerFlux(Layer="top")
    except Exception as e:                                   # noqa: BLE001
        return str(e), ["zero"]
    return None, ["zero"]


def case_left_handed_lattice_ok():
    """A negatively oriented basis is a valid cell and must not be refused."""
    S = s4compat.New(Lattice=((1, 0), (0, -1)), NumBasis=25)
    S.SetMaterial(Name="bg", Epsilon=1.0)
    S.SetMaterial(Name="e9", Epsilon=9.0)
    S.AddLayer(Name="top", Thickness=0, Material="bg")
    S.AddLayer(Name="pat", Thickness=0.35, Material="bg")
    S.SetRegionCircle(Layer="pat", Material="e9", Center=(0, 0), Radius=0.2)
    S.AddLayer(Name="bot", Thickness=0, Material="bg")
    S.SetExcitationPlanewave(IncidenceAngles=(0, 0), sAmplitude=1, pAmplitude=0)
    S.SetFrequency(1 / 1.6)
    try:
        S.GetPowerFlux(Layer="top")
    except Exception as e:                                   # noqa: BLE001
        return str(e), None
    return None, None


def case_ellipse_is_a_circle():
    """An ellipse with equal half-widths is a circle, at every angle.

    Two exact identities, neither of which depends on knowing the right answer:
    an ellipse with equal half-widths is a circle and rotating it changes
    nothing, and an ellipse rotated by 90 degrees is the same shape with its
    half-widths swapped.

    The ellipse branch of shape_get_tangent_cross_segment -- which the normal
    vector field is built from, so it runs for every polarization basis -- had
    three faults at once: the segment direction was rotated into the ellipse
    frame using only its x component, the intersection loop indexed by the loop
    count rather than the loop variable (reading isect[4] out of a four-element
    array when there were two intersections), and the dot product used the
    unscaled component.  Upstream still has all three.

    Before the fix the ellipse misses the circle by 4e-2 and is not itself
    rotation-invariant -- 0 degrees sits 1.6e-2 from the other angles, because
    sa vanishes there and only one of the two wrong components costs anything.
    Rotation is not being discarded: an ellipse of unequal half-widths responds
    to it either way.
    """
    def transmission(place, basis, nb=100):
        S = s4compat.New(Lattice=((1, 0), (0, 1)), NumBasis=nb)
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

    centre = (0.05, -0.03)                      # off-centre, so rotation bites
    circle_at = lambda r: (lambda S: S.SetRegionCircle(                 # noqa: E731
        Layer="pat", Material="rod", Center=centre, Radius=r))
    ellipse_at = lambda a, b, ang: (lambda S: S.SetRegionEllipse(       # noqa: E731
        Layer="pat", Material="rod", Center=centre, Angle=ang, Halfwidths=(a, b)))

    for basis in ("Normal", "Default"):
        want = transmission(circle_at(0.30), basis)
        if want != want:                        # NaN
            return "circle gives NaN under %s" % basis, None
        for angle in (0, 17, 90):
            got = transmission(ellipse_at(0.30, 0.30, angle), basis)
            if got != got:
                return "ellipse at %d degrees gives NaN under %s" % (angle, basis), None
            if abs(got - want) > 1e-11:
                return ("ellipse(0.3,0.3) at %d degrees differs from the circle by "
                        "%.2e under %s" % (angle, abs(got - want), basis)), None
        a, b = 0.34, 0.18
        one = transmission(ellipse_at(a, b, 25), basis)
        two = transmission(ellipse_at(b, a, 115), basis)
        if abs(one - two) > 1e-11:
            return ("(%.2f,%.2f) at 25 and (%.2f,%.2f) at 115 differ by %.2e under %s"
                    % (a, b, b, a, abs(one - two), basis)), None
    return None, None


CASES = [
    ("crossing regions",            case_cross,             "refuse"),
    ("region meets its own image",  case_self_image,        "refuse"),
    ("region meets another image",  case_other_image,       "refuse"),
    ("unknown material",            case_unknown_material,  "refuse"),
    ("1D overlapping intervals",    case_1d_overlap,        "refuse"),
    ("1D nested intervals",         case_1d_nested,         "refuse"),
    ("1D pair across the boundary", case_1d_across_boundary, "refuse"),
    ("1D wider than the cell",      case_1d_wider_than_cell, "refuse"),
    ("nested regions",              case_nested_ok,         "accept"),
    ("region spanning the cell",    case_spanning_ok,       "accept"),
    ("1D disjoint intervals",       case_1d_disjoint_ok,    "accept"),
    ("1D region filling the cell",  case_1d_full_cell_ok,   "accept"),
    ("clone of a patterned layer",  case_clone,             "accept"),
    ("both keyword spellings",      case_both_keyword_spellings, "accept"),
    ("long non-ASCII layer name",   case_long_name_message, "refuse"),
    ("SetVerbosity",                case_set_verbosity,     "accept"),
    ("version and GetOptions",      case_introspection,     "accept"),
    ("docstrings name themselves",  case_docstrings,        "accept"),
    ("plane interface = Fresnel",   case_two_layer_interface, "accept"),
    ("duplicate layer name",        case_duplicate_layer_name, "refuse"),
    ("degenerate lattice",          case_degenerate_lattice, "refuse"),
    ("1D period of zero",           case_zero_period_1d,    "refuse"),
    ("SetLayer redefines",          case_redefine_layer_ok, "accept"),
    ("left-handed lattice",         case_left_handed_lattice_ok, "accept"),
    ("ellipse(r,r) is a circle",    case_ellipse_is_a_circle, "accept"),
]


def main():
    print("S4 build :", s4compat.where())
    print("keywords :", s4compat.STYLE)
    print()
    failures = 0
    for label, fn, expect in CASES:
        try:
            msg, want = fn()
        except Exception as e:                               # noqa: BLE001
            print("  FAIL  %-28s case itself raised %s: %s"
                  % (label, type(e).__name__, e))
            failures += 1
            continue

        if "accept" == expect:
            if msg is None:
                print("  ok    %-28s accepted" % label)
            else:
                print("  FAIL  %-28s should have been accepted: %s" % (label, msg))
                failures += 1
            continue

        if msg is None:
            print("  FAIL  %-28s accepted, but is not representable" % label)
            failures += 1
            continue
        missing = [w for w in (want or []) if w not in msg]
        if missing:
            print("  FAIL  %-28s refused, but the message omits %s"
                  % (label, ", ".join(repr(m) for m in missing)))
            print("        %s" % msg)
            failures += 1
        else:
            print("  ok    %-28s refused: %s" % (label, msg))

    print()
    print("%d/%d cases behave as intended" % (len(CASES) - failures, len(CASES)))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.path.insert(0, __file__.rsplit("/", 1)[0])
    sys.exit(main())
