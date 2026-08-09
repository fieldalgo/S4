# Validation report: the Python binding port

Whether the Python binding port on `modernize/master` changed the physics, and
what the checking turned up on the way. There is one criterion: **put the same
input through the 2016 oracle and through the new build, and compare the numbers
digit by digit.** Quantities the oracle cannot produce are pinned instead by
self-consistency -- energy conservation, lattice periodicity, symmetry, the sum
over diffraction orders.

- Baseline: benvial `ad4f128`, the build installed in the conda env `s4` (2016 physics)
- Subject: `modernize/master` plus the branches below
- Where the mechanisms are: each defect's cause and its smallest reproduction
  live in `UPSTREAM.md`, written to be sent to victorliu/S4. This report records
  what was measured, what it establishes, and what it does not.
- Harness: `validation/` (70 probes and 24 behaviour cases; the port-verification
  table in §2 was taken at 62, before the Priority 4 and later probes were added)

---

## 0. Conclusion first

**The port itself is sound.** Every probe the oracle can be compared against
agrees to a relative error of 1e-12 or better. Four disagreed; all four are
explained, and in three of them the new build is the correct one.

**Fourteen defects unrelated to the port were found.** All of them are in
upstream master's C core or its binding, eight produce **silently wrong answers**
and one is a **SIGSEGV**.

The heaviest is the simplest. **A structure of exactly two layers -- a plane
interface, the only case with a closed-form answer -- comes back with T = 0,
R = 0 and R+T = 0**, silently. Three layers or more are unaffected, which is how
it survived, and why the 62-probe sweep in §2 missed it: every probe built at
least three.

Next after it is the FMM grid. **With no options touched at all**, a large enough
basis (`NumBasis=2453`, say) lets uninitialised memory into the permittivity
coefficients. The answer can change from run to run; `R+T = 4.5e+133` was
observed once.

| Defect | Symptom | Where | Branch |
|---|---|---|---|
| **two-layer structure** | a plane interface returns T = 0, R = 0, R+T = 0 | `S4/rcwa.cpp` | `fix/two-layer-and-refusals` |
| ellipse tangent code | rotation dropped, loop indexed by its count, read past the array; every basis affected | `S4/pattern/pattern.c` | `fix/ellipse-tangent-constraint` |
| `Clone()` | using the clone segfaults, and the source is damaged too | `S4/S4.cpp` (upstream, from the v2 API transition) | `fix/simulation-clone` |
| rotation-centre typo | nesting inside a rotated, off-diagonal-centre parent disappears | `S4/pattern/pattern.c:159` | `fix/pattern-rotation-centre` |
| crossing regions | one region is dropped whole without warning (2D); the shared part is double-counted (1D) | `S4/pattern/pattern.c`, `S4/S4.cpp` | `fix/detect-crossing-regions` |
| region meeting a periodic image | shared area counted once per copy; mean permittivity exceeds what the materials allow | `S4/S4.cpp` | `fix/periodic-image-overlap` |
| **discretisation grid** | wrong values on an odd `ngrid`; a read past the buffer at `res=2` | `S4/fmm/*.cpp` | `fix/fmm-discretized-grid` |
| **`LatticeTruncation` ignored** | asking for `Parallelogramic` leaves the basis circular. No error. The count sometimes matches by chance | `S4/main_python.c` + when G is selected in `S4.cpp` | `fix/lattice-truncation` |
| `SolveInParallel` | an empty stub | `S4/main_python.c` | `fix/module-entry-points` |
| `NewSpectrumSampler` | malformed format string → `SystemError` | `S4/main_python.c` | `fix/module-entry-points` |
| `SetVerbosity` | inverted parse check: unreachable from Python, `SystemError` on every call | `S4/main_python.c` | `fix/setverbosity` |
| duplicate layer name | both layers change the answer, only the first can be addressed | `S4/S4.cpp` | `fix/two-layer-and-refusals` |
| degenerate lattice | parallel basis vectors accepted; NaN results, LAPACK complaint on stderr | `S4/S4.cpp` | `fix/two-layer-and-refusals` |
| `GetEpsilon` return protocol | returns a value with an exception set → overwritten by `SystemError` | `S4/main_python.c` | `fix/detect-crossing-regions` (first commit) |
| rectangle normal sign `-.1` | **no observable effect** (normalisation absorbs it) | `S4/pattern/pattern.c:248,251` | `fix/pattern-rectangle-normal-sign` |

The last two change no result; they are listed because they are real and cheap
to fix, not because they cost anything today.

---

## 1. The harness

Everything lives in `validation/` inside the S4 repository. Nothing was placed
in any sibling project's tree.

| File | Role |
|---|---|
| `probe.py` | **actually calls** every binding entry point and emits numbers. Asserts nothing |
| `compare.py` | runs the two builds in separate processes and compares digit by digit. A probe that segfaults is recorded, not fatal |
| `diagnostics.py` | asserts behaviour rather than numbers: what must be refused, what must keep working, and what the message says |
| `s4compat.py` | absorbs the one API difference between the builds (`Layer=` ↔ `S4_Layer=`) |
| `build.sh` | isolated build and install, in two flavours. **Guards the shared env by hash** |
| `normals_probe.c` | links `shape_get_normal` directly out of libS4.a, which Python cannot reach |

```bash
sh validation/build.sh /tmp/s4-mine                  # library plus isolated install
python validation/compare.py --a "" --b /tmp/s4-mine \
       --label-a oracle --label-b mine               # --a "" = the oracle installed in the env
python validation/compare.py --a /tmp/s4-a --b /tmp/s4-b  region_nested dc_crossing_regions
```

`build.sh` never runs `make S4_pyext`: that target's last line is
`pip3 install --upgrade ./`, which overwrites whatever environment is active. It
compares the sha256 of the
guarded install (`S4_GUARD`, an extension in a shared environment that must not
be touched) before and after, and fails the build if it moved. Throughout this work it stayed
`81c4a0946e0382d602dd9db3a6a7ac2d5aa714b7c73e29682997f5b1297d6a51`.

### How the geometry is measured: the DC Fourier coefficient

`GetEpsilon` returns the band-limited permittivity S4 actually solves. Averaged
on a grid past Nyquist, that mean **is exactly the DC coefficient** -- the
truncation error is zero by construction, not small. So the quantity separates
"the geometry went in wrong" from "it has not converged yet". Every geometric
verdict below is this quantity.

---

## 2. Priority 1 -- verifying the port

`validation/probe.py`, 62 probes, oracle vs `validate/binding-sweep` (`d03ca3e`).

```
same: 55   differs: 4   failed: 3
```

### 2.1 What agreed (55/62)

The main entries:

| Area | What was checked |
|---|---|
| reference convergence series | NumBasis 200/400/600/800 → `0.842114450429` / `0.840432203521` / `0.840907273718` / `0.840475065200`. Relative error 6e-14 to 6.8e-13 |
| four bases | no decomposition / Default / Normal / Jones all agree. Jones's `R+T=1.0000016813` excursion matches the oracle to the last digit too |
| four region types | circle, ellipse, rectangle, polygon. Rectangle = the same polygon (2.6e-14); a 90° rotation = swapped half-widths (2.0e-13) |
| nested / disjoint / ring / four siblings | DC coefficient within 1e-15 of the analytic value |
| **layer copies** | `AddLayerCopy` = an explicitly duplicated layer (`diff` 1e-14). A plain layer accepts patterning, a copy refuses it -- the predicate `b6fe76b` fixed |
| **diffraction orders** | `GetBasisSet` returns 113, origin included, centrosymmetric. At λ=0.7, per-order transmission plus `Σ orders = total` to 1e-16 |
| **off-normal / conical incidence** | θ=10/30/55°, s and p separately; (θ,φ)=(30,30),(30,45),(45,20). On a C4 structure φ and φ+90° agree |
| **anisotropy** | diagonal tensor = scalar; uniaxial tensor s≠p; **off-diagonal complex tensor `T=0.866227334716`** (so the zgeev fix really works); gyrotropic (Hermitian) `R+T=1`; a tensor material inside a pattern |
| **field output** | `GetFields` matches the analytic plane wave in a uniform slab; lattice translation periodicity 2e-15; C4 rotational symmetry; `GetFieldsOnGrid` ↔ `GetFields` |
| `GetAmplitudes` | the incident layer's order-0 forward amplitude is exactly the input, the rest of the forward amplitudes are 0 |
| `GetStressTensorIntegral` / `GetLayerVolumeIntegral` (U,E,H,e) / `GetLayerZIntegral` | all agree |
| `SetExcitationExterior` | order-0 s excitation = plane-wave excitation (`diff` 3e-14). Non-specular orders, y polarisation and mixtures agree too |
| complex frequency | `SetFrequency(complex(1/1.6, -1e-3))` agrees, and the imaginary part really takes effect (T moves by 0.0040) |
| 1D lattice | NumBasis 11/31/51/101 convergence, 1D ↔ 2D stripe agreement, order decomposition at λ=0.7, θ=0/20/40°, 1D fields |
| the rest | `GetSMatrixDeterminant` (outside `Clone`), `GetReciprocalLattice`, hexagonal lattice, `RemoveLayerRegions`, redefinition via `SetLayer`, absorbing materials, file output |

### 2.2 The four that disagreed -- all accounted for

**(a) `region_circle_vs_ellipse` -- the new build is right.**

| | oracle | new build |
|---|---|---|
| circle R=0.30 | 0.846175819327 | 0.846175819327 |
| ellipse (0.30, 0.30) | **0.826151184323** | **0.846175819327** |

The oracle's ellipse is wrong. This is the difference `f020606` (the ellipse
tangent-constraint fix) intended, and in the new build circle and ellipse agree
to 1.6e-14.

Revisited later with an identity rather than a single number, which says more.
An ellipse with equal half-widths **is** a circle, and rotating a circle changes
nothing; that needs no reference implementation to check. Off-centre, with
`PolarizationDecomposition` on:

| shape | oracle | fixed |
|---|---|---|
| circle, R = 0.30 | 0.807627216212 | 0.807627216212 |
| ellipse (0.30, 0.30) at 0 deg | 0.782190513141 | 0.807627216212 |
| ellipse (0.30, 0.30) at 25, 45, 70, 90 deg | 0.765815666420 (all four) | 0.807627216212 |

Two failures at once on the oracle: writing the circle as an ellipse moves the
answer by 4.2e-2, and **the ellipse is not rotation-invariant** although a
circle must be -- 0 degrees is 1.6e-2 away from every other angle. Under `Jones`
and the default basis two of the cases return `NaN`. The fixed build holds all
five angles to 6.4e-14 and matches the circle.

Zero is special because `sa` vanishes there, so the term missing from `u[0]`
costs nothing and only `u[1]` is wrong; at any other angle both components are.
Why the remaining angles agree with one another was not traced.

Rotation itself is not being discarded -- an ellipse of unequal half-widths
moves T by 0.15 across these angles on both builds. The second identity says
the same thing from the other side.

**Isolated properly.** The comparison above is against the 2016 oracle, which is
a different codebase, so it carries every 2016-to-2018 difference with it. To
separate the ellipse fix alone, `f020606` was reverted on top of an otherwise
identical `modernize/master` build. Of 70 probes, 2 change --
`region_circle_vs_ellipse` and `discretization_grid_parity`, which builds an
ellipse -- and everything else is bit-identical.

| condition | affected |
|---|---|
| `PolarizationDecomposition` off | no -- both builds identical, ellipses included |
| circle, rectangle, polygon | no -- bit-identical |
| ellipse, equal half-widths | yes, 2.5e-2 to 4.2e-2 |
| ellipse, unequal half-widths | yes, 1.5e-2 at 0 degrees and 4.0e-2 at 90 |
| either, at 0 degrees | **non-deterministic** -- `NaN`, 0.49, 0.84 and 1 across eight runs of one command |

The last row is fault (2) showing itself: `isect[4]` is one element past a
four-element array, so the answer depends on what the stack happened to hold.
LAPACK prints `ZGEEV parameter number 4 had an illegal value` to stderr each
time. The fixed build returns the same numbers in all eight runs.

Away from 0 degrees the three polarization bases agree with each other to 3e-8
on the unfixed build, where after the fix they differ by up to 1.8e-2 as
different factorisations should -- the decomposition was being built from a
field that no longer told them apart.

That the tangent code is reached only through `Pattern_GenerateFlowField`, which
only the `PolBasis*` formulations call, predicts the first row; measuring it
confirms the reasoning rather than assuming it.

**How far the fix is verified.** The circle identity is exact and the circle
takes a different branch of the same function, so it is an independent route.
Approximating the ellipse by an N-gon -- a third route, through
`intersection_polygon_segment` -- converges onto the fixed value for that case
(5.5e-3 at 24 sides to 3.9e-5 at 384) and never approaches the unfixed one.

For unequal half-widths there is no identity, and the polygon route cannot
substitute for one: under the normal-vector basis it plateaus 1.7e-3 and 6.2e-3
away and stops moving between 384 and 768 sides. With the decomposition off the
same refinement converges cleanly, a factor of four per doubling, which places
the residue in the normal field rather than the geometry -- a faceted polygon
keeps its corners however fine it gets. So the fix is confirmed exactly where an
identity exists, and elsewhere only shown to be much closer to an independent
construction than what it replaced.

The cause is three faults in one branch of `shape_get_tangent_cross_segment`,
which the normal vector field is built from and which therefore runs for every
polarization basis: the segment direction is rotated using only its x component,
the intersection loop subscripts by the loop *count* rather than the loop
variable -- reading `isect[4]` out of a four-element array when there are two
intersections -- and the dot product uses the unscaled component. Upstream still
has all three; UPSTREAM.md 14 writes it up.

The companion identity, an ellipse rotated 90 degrees against the same ellipse
with swapped half-widths, holds to 8e-15 even before the fix. The geometry was
never the problem; the direction handling was. Afterwards both identities hold
to 1.0e-13 across all three bases, which `diagnostics.py` asserts.

**(b) `fields_on_grid` -- the array index convention is transposed. The values are the same.**

Both builds agree with `GetFields` to 1e-15, but **the axis order is reversed**.

| | the point `E[i][j]` refers to |
|---|---|
| oracle (2016) | `(x, y) = (j/n, i/n)` -- the first index is y |
| new build (2018) | `(x, y) = (i/n, j/n)` -- the first index is x |

**This is where a harness silently transposes when it is moved.** It is an
array convention rather than a physical quantity, so no amount of energy
conservation will catch it. The probe measures both readings and records which
one holds.

> Note: this grid is an inverse FFT of the mode amplitudes. If the grid cannot
> hold the highest harmonic in the basis, disagreeing with `GetFields` is
> correct behaviour. NumBasis 80 (|m|≤5) needs n≥11. Measured at n=8 the two
> builds are wrong in different ways, and that is aliasing, not a defect.

**(c) `amplitudes_incident` -- length 240 vs 226. The new build is self-consistent.**

`GetAmplitudes` returns an array of length 240 (=2×120) from the oracle and 226
(=2×113) from the new build -- while both builds' `GetBasisSet()` returns
**113**. So the oracle sizes the buffer from the requested 120, fills only the
113 that exist, and reports the buffer's length. The values themselves (order 0
= −1, the rest of the forward amplitudes 0, backward maximum 0.392204246913)
agree exactly.

**(d) `output_files` -- `fld.E` 1796 vs 1797 bytes.** The transposition in (b)
changes the `%.14g` digit count by one character. The other four files match to
the byte.

### 2.3 The three failures -- `Clone()`

`clone_matches` / `clone_full_structure` / `clone_is_independent` **SIGSEGV** on
the new build and are fine on the oracle. §4.

### 2.4 So, the port?

**The pointer → integer-handle port introduced no physical difference.** Every
path the oracle can answer on agrees to the last digit. The functions the
briefing worried about -- the ones that would link but break when called --
were all actually called and all work: `GetFields`, `GetFieldsOnGrid`,
`GetAmplitudes`, `GetPowerFluxByOrder`, `GetStressTensorIntegral`,
`GetLayerVolumeIntegral`, `GetLayerZIntegral`, `GetBasisSet`,
`SetExcitationExterior`.

Many of them still use the **old pointer API**
(`Simulation_GetLayerByName(S, name, NULL)`). They link and run because master
never removed that internal API. The port moved only the seven functions that
had disappeared over to v2, and that judgement was right.

---

## 3. Priority 2 -- two typos

### 3.1 The rotation centre (`shape_contains_point`) -- **a real physics bug**

```c
x[1] = (x_[0] - s->center[1]) *-sa + (x_[1] - s->center[1]) * ca;
                        ^^^^ should be center[0]
```

The query point is displaced by `sin(angle) * (center[1] - center[0])`, a
constant independent of the query point. It is 0 without rotation, and 0 on the
diagonal (`cx == cy`) as well. So a rotation test centred on the origin, or an
off-centre test without rotation, can never catch it -- and neither can a layer
with a **single** region, because the containment tree does not run at all.

The predicate is used by `pattern_get_containment_tree`. A displaced query
answers about the wrong point, so a region properly inside a rotated,
off-diagonal-centre parent **is recorded as sitting on the background**, and its
Fourier contribution is computed against the wrong reference permittivity. This
is not an accuracy question: **a different structure is solved.**

Measured (`dc_rotated_offcentre_nesting`): an eps-9 rectangle centred at
(0.15, −0.15), rotated 45°, half-widths (0.40, 0.08), holding an eps-4 circle of
radius 0.04. The displacement is 0.212 against a half-width of 0.08.

| | mean eps |
|---|---|
| exact | **1.998867258770** |
| S4 before | **2.039079644740** ← matches, to 3.1e-15, the value for a circle sitting on the background |
| S4 after | 1.998867258770 (error 3.1e-15) |

The control `dc_rotated_offcentre_control` is the same geometry with `Angle=0`.
The displacement is identically zero there, and it is accurate to 4.4e-16 both
before and after -- so nothing but the rotation can explain the difference.

**Every other probe is identical before and after**: 58 agree digit for digit,
and 3 are `Clone` and die the same way on both sides (§4, unrelated to this
branch).

Branch `fix/pattern-rotation-centre`, commit `cd9f4a9`.

### 3.2 The rectangle normal sign `-.1` -- **no observable effect**

```c
double sgn = (rx > 0. ? 1. : -.1);   /* two places */
```

An obvious typo that **changes no result.** `shape_get_normal` normalises just
before returning, so a vector built with length 0.1 comes back as the unit
vector in the same direction. All that survives is rounding.

Read straight out of libS4.a with `validation/normals_probe.c` (a rectangle at
three angles × four faces = 12 components):

| | before | after |
|---|---|---|
| rotated 57°, −x face | `(-0.54463903501502708, -0.83867056794542394)`, \|n\| = 0.99999999999999989 | `(-0.5446390350150272, -0.83867056794542405)`, \|n\| = 1 |

**Two of the twelve components move, by 1 ulp each.** The Python-level probes
are **bit-identical**: 59 agree and the remaining 3 are `Clone`, dying the same
way on both sides.

So this is not a bug fix but the removal of a trap. The only current consumer is
the Kottke subpixel branch in the FMM, which renormalises anyway. A future caller
that trusts the documented unit-length contract and does not renormalise would
get a vector ten times too short on two of the four faces.

Branch `fix/pattern-rectangle-normal-sign`, commit `470ef1b`.

---

## 4. Turned up by Priority 1 -- `Clone()` segfaults

On `modernize/master`, using a clone kills the process. The clone pretends to
succeed and the first lookup of a layer or material by name blows up. The 2016
oracle is fine. **This is not the port's doing but a defect in upstream master's
C core** -- the 2016 version wrote into array slots directly, and the move to the
v2 append API broke it.

There are four causes, layered.

1. `memcpy(T, S, sizeof *S)` copies **the counts as well**, and then fills the
   freshly allocated arrays with `id = -1` (append). Append goes after the
   current count, so the real entries land in `[n, 2n)` and `[0, n)` stays as
   malloc gave it. A name lookup scans from 0 and calls `strcmp` on garbage
   pointers.
2. That setter is handed `M->type` (the internal 0/1 encoding) where it expects
   `S4_MATERIAL_TYPE_*` (which starts at 2). It always falls to `default`, which
   **removes the material just added and copies no epsilon at all.** Even with
   cause 1 fixed, the clone would have had no materials.
3. `S4_Simulation_SetLayer` calls the target's `Simulation_DestroySolution`.
   `T->solution` still points at `S->solution` from the memcpy, so adding the
   first layer to the clone **frees the source's solution** while the source goes
   on pointing at it.
4. `G`, `kx`, the vector-field dump prefix and the polygon vertex arrays stay
   shared, so destroying both is a double free.

The fix zeroes the counts and the two cache pointers **before** the loops, deep-
copies `G`, `kx`, the prefix and the polygon vertices, and passes the complex
type constants (`eps` is a union overlaying `s[2]` and `abcde[10]`, so the
complex form contains the real form exactly). Materials and layers are appended
to empty arrays in the source's order, so ids are preserved and the `copy`,
`material` and `exc.layer` references stay valid.

Verification: three probes, SIGSEGV before, and after:

| | |
|---|---|
| `T_orig` vs `T_clone` | `diff` = **0** (identical) |
| `source_damage` (the source re-measured after cloning) | **0** |
| editing only the clone, then the source | `source_damage` = 0, the clone moved by 0.148 |
| polygon + tensor + layer copy + conical excitation | `T` and `eps` both `diff` = 0 |

**The other 59 probes are unchanged** (differs 0).

Branch `fix/simulation-clone`, commit `6bbac7f` (branch head `6bbac7f`).

---

## 5. Priority 3 -- crossing regions

### 5.1 Confirmation

The briefing's numbers reproduce. Two squares of half-width 0.17 (eps 4 and
eps 9) at x = ∓0.08.

| | mean eps |
|---|---|
| exact (the later region wins) | **2.0880** |
| S4 | **1.9248** |
| the value when eps-4 is forced to be eps-9's parent | 1.9248 (error 2.0e-15) |

So **the eps-4 square's exposed area vanishes entirely.** Without a warning.
The legitimate cases are exact: single 1.3468, disjoint 2.2716, nested 5.5420,
all within 1e-15.

### 5.2 Why it was not caught, and what replaced it

`pattern_get_containment_tree` contains the code meant to prevent this, and two
things stopped it: a pair loop that advances the wrong variable, and a
`shapes_intersect` whose every branch returns 0 (UPSTREAM.md 3).

The replacement is two complementary tests, either of which is enough to refuse.
`boundaries_cross` is exact -- transversal segment crossing for polygons, the
quadratic in the conic's normalised frame for polygon-conic, a sign-change sweep
for conic-conic. `shapes_overlap_partially` samples both regions and refuses when
each has interior points inside and outside the other, which is what catches two
identical axis-aligned rectangles offset along one axis: the commonest form of
the mistake, whose boundaries meet only along collinear edges and so never cross.

Both are **biased towards acceptance**. Every refusal has a witness -- a
transversal crossing, a sign change, or a point inside one and outside the other
-- so nothing S4 accepts today is newly refused unless it really does overlap.
That bias is what §8.1a later showed to be imperfect, and what the tolerance fix
restored.

### 5.4 1D was worse -- **nesting is silently wrong too**

Not in the briefing. A 1D region is a rectangle with half-width[1] = 0, so its
**area is zero** and `shape_contains_point` accepts no point at all. The
containment tree therefore parents every 1D region to the background and the
transform adds each interval's contribution as it stands.

| 1D case | exact | S4 before |
|---|---|---|
| single (half-width 0.25) | 2.5 | 2.5 |
| disjoint | 3.2 | 3.2 |
| abutting (shared endpoint) | 3.2 | 3.2 |
| **nested** (0.15 inside 0.40) | **5.9** | **8.3** -- wrong |
| **partial overlap** | **4.8** | **5.4** -- wrong |

There is no notion of nesting in 1D at all. The 2D distinction -- nesting is
fine, crossing is unrepresentable -- does not hold. So **overlapping 1D
intervals are refused whether the overlap is partial or total** (decided in
`S4.cpp`, which knows the lattice). Intervals that share only an endpoint are
allowed; that is how a nested structure has to be written in 1D.

### 5.5 Getting the error to the caller

The containment tree encodes which region is at fault in the magnitude of its
return value. Handing that back verbatim collides with `Simulation_InitSolution`'s
own small error codes -- **a layer with two regions reported "A memory allocation
error occurred".** It is folded into code 17 (invalid region) and 18 (overlap),
with entries added to both frontends' message tables.

Along the way, `GetEpsilon` turned out to **return a value with an exception
set.** CPython detects this and overwrites the original error with
`SystemError: ... returned a result with an exception set`. Every other entry
point in the file returns NULL. Separate commit `a08da3c`.

### 5.6 Thin overlaps -- strengthened afterwards (`ae1428e`)

The first test accepted two axis-aligned rectangles overlapping by 0.01 or less.
A rectangle offset along one axis has collinear top and bottom edges, where
`boundaries_cross` deliberately declines to answer so that touching shapes are
not refused; the interior-point test that fills the gap was scattering its 32x32
grid over each shape's whole bounding box, whose spacing steps over a thin band.

Moving only the "inside both" search to the **intersection** of the two bounding
boxes makes the grid spacing shrink with the overlap. It now catches down to
1e-8, where before 1e-2 passed.

**Still biased towards acceptance**, and the probe is built to show it: every
overlap case is measured alongside the same shapes moved apart by the same
distance. All six gap cases down to 1e-5 pass before and after -- a test that
refuses everything scores full marks on the overlap half alone.

```
overlaps caught  6/10 -> 10/10
gaps preserved   6/6  -> 6/6
```

### 5.7 Verification

`crossing_detection` -- **all ten unordered pairs** of the four region types,
once overlapping and once legitimate, plus triple nesting, four siblings,
external and internal tangency, and a rotated off-centre parent: 21 patterns.

| | before | after |
|---|---|---|
| the 11 overlapping | all accepted (silently wrong) | all `RuntimeError` |
| the 10 legitimate | all accepted | all accepted (unchanged) |

`crossing_detection_1d` -- the table in §5.4. The two overlapping cases refused,
the three legitimate ones preserved to 1e-15.

**This branch has to sit on top of `fix/pattern-rotation-centre`.** The
interior-point test uses `shape_contains_point`, so with the rotation-centre typo
still present it refuses **legitimate** nested pairs that have a rotated,
off-diagonal-centre parent. Measured, not assumed.

Branch `fix/detect-crossing-regions`, commits `a08da3c` + `43e266d` + `ae1428e`
(branch head `35a37d7`).

### 5.8 Except the message was being thrown away (`fix/error-detail`)

§5.5 says the layer name was put into the message. **It did not reach the
caller.**

The Python binding never calls `S4_Simulation_SetMessageHandler`, so `S->msg` is
NULL and the carefully built `"Layer 'my_pattern' has two regions..."` is
discarded whole. What was left is the one generic sentence from the error-code
table.

```
GetPowerFlux: A layer has two regions whose boundaries cross; regions must be
either nested or disjoint
```

True, and useless on a structure with forty layers. This hole was of my own
making.

**Three things fixed.**

**(1) Carrying the message to the exception.** `S4Sim` installs a handler that
keeps the most recent error message, and `HandleSolutionErrorCodeDetail`
prepends it. The table entries were rewritten to carry the **remedy** rather
than restate the diagnosis, since the diagnosis now arrives ahead of them.

**(2) Which region.** The containment tree already encoded the offending index in
its return value (§5.5), and the periodic-overlap check returned one too. Both
were being discarded and replaced with "a layer has ...". Both checks run
**before the sort by area**, so the index is the order the caller added the
regions in -- the only numbering the caller ever sees.

The periodic-overlap message had to name **both** members of the pair. A small
region near the edge landing inside a large region's copy is commoner than a
region meeting its own repeats, and the message claimed "its own" in both cases.

```
GetPowerFlux: Region 3 of layer 'grating' crosses another region.
  Regions must be either nested or disjoint; they are numbered from 1 in the
  order they were added to the layer

GetPowerFlux: Region 1 of layer 'grating' overlaps the periodic image of
  region 2. ...
```

**(3) The 1D check could not see the period.** The 1D interval check added in
§5.4 subtracted centre coordinates directly. It saw overlaps written inside the
cell and **not overlaps that meet across the boundary.**

| 1D structure | before | after |
|---|---|---|
| centres −0.45 and +0.45, half-width 0.1 (meeting at the boundary) | accepted (silently wrong) | refused |
| a single half-width 0.6 (wider than the period) | accepted (silently wrong) | refused |
| a single half-width 0.5 (exactly filling the cell, touching only) | accepted | accepted |

The 2D periodic check cannot catch these: a 1D region has zero height and
intersects nothing. The separation is now reduced modulo the period, and a region
wider than the period is refused outright.

The 1D refusal uses **a new code 20** rather than 18, because 18's remedy is
"nested or disjoint" and nesting is not an option in 1D at all -- the 1D
transform **adds** each interval against the layer background and never subtracts
the enclosing one. An eps-4 half-width 0.1 inside an eps-9 half-width 0.3 reads a
mean permittivity of **6.4**, the naive superposition `1 + 8·0.6 + 3·0.2`, where
the truth is 4.8.

Branch `fix/error-detail`, commits `899e0a9` + `3d4a036` + `685ed8c`.

### 5.9 A long or non-ASCII layer name ate the complaint (`fix/message-truncation`)

The messages are built in a 256-byte buffer that the layer name is free to fill.
When it does, `snprintf` cuts the sentence instead of the name -- and the
sentence is where the complaint is. A layer named with eighty Hangul syllables
produced the name followed by the generic remedy, with `crosses another region`
gone: the reader is told the numbering convention and never told which region.

The cut is also byte-wise, so for any name that is not ASCII it lands inside a
multi-byte character often enough to matter, and Python renders the remains as
U+FFFD. Both failures came out of the same call.

The name is elided instead, to 96 bytes and on a character boundary, with a
marker so it is visible that it happened. Stepping back off continuation bytes
(0x80-0xBF) reaches a character start for any input, well formed or not.

| layer name | before | after |
|---|---|---|
| ASCII, short | complete | complete |
| 80+ non-ASCII characters | U+FFFD, diagnosis truncated away | elided with `...`, diagnosis intact |

Branch `fix/message-truncation`, commit `dd767a1`.

### 5.10 Verification of §5.8 and §5.9: `validation/diagnostics.py`

Every other script here compares **numbers** between two builds. These defects do
not appear as a wrong number in a converged series -- they appear as a plausible
number for a structure the caller did not describe. So they need a test that
compares **behaviour**.

Twenty-five cases: nine structures that cannot be represented and must be
refused, one lookup failure, and fifteen behaviours that must keep working --
including the plane interface of §9, which is asserted against the closed form
rather than against another build. Each
refusal also asserts what the message says. The accepted cases are what stops the
checks from being tightened into uselessness, and they include the two shapes
closest to the line -- a nested pair in 2D and a region exactly as wide as the
cell.

```
new build    25/25
2016 oracle   9/25
```

The oracle's fifteen failures are the inputs it accepts and then solves as
something else, plus the entry points it cannot offer at all. The 70 probes are
**entirely unchanged** by the diagnostics work
(`same: 70  differs: 0  failed: 0`).

---

## 6. Priority 4 -- `res=4` non-determinism: **diagnosed and fixed**

### 6.1 Reproduction: non-deterministic on the oracle, mostly deterministic on the new build

The same input (NumBasis 400, `DiscretizedEpsilon`,
`DiscretizationResolution=4`, Normal basis) run 20 times on each.

| build | result |
|---|---|
| 2016 oracle | **three different answers** -- 12× `T=0`, 5× `T=0.863944858`, 3× `NaN` |
| `validate/binding-sweep` | 20/20 `T=0.863944858` |

On the oracle LAPACK complains directly:

```
** On entry to  ZGEEV parameter number  4 had an illegal value
Layer eigensystem returned info = -4
```

So the oracle's non-determinism is **already-corrupted values reaching the
eigensolver**. After the rwork fix (`21302a5`) that error disappears and it looks
deterministic.

**But "deterministic" is an illusion.** The defect in §6.3 reads uninitialised
heap, so the answer depends on what was allocated earlier in the same process.
During this work, `nb=120, res=2` produced `R+T = 4.5e+133` once, where the same
command usually gives `1.00006414`. Same character as the divergence recorded in
the briefing (1.06e+111), and **it still happens on the current build.** It just
does not reproduce on demand -- which is the nature of the defect.

### 6.2 The two causes

Both are in UPSTREAM.md 5 and 6. In brief: the half-grid rotation
`si = ii >= n/2 ? ii - n/2 : ii + n/2` is not a bijection when `n` is odd, so a
row and a column of the grid keep whatever `S4_malloc` returned; and the
coefficient read wraps only the negative side of `f = G_a - G_b`, so it indexes
outside the buffer whenever `ngrid == 2*Gmax`, which `res=2` reaches exactly.

The damage is not confined to high frequencies. NumBasis 400, res=4 (grid 45):

| | DC coefficient error | worst coefficient error (location) |
|---|---|---|
| before | **4.395e-02** | 2.776e-01 at f=(0,2) -- got -0.229, exact **+0.042** |
| after | **7.99e-15** | 3.000e-01 at f=(1,1) -- the same place as an even grid |
| res=5 (grid 60, even) | 1.33e-14 | 1.129e-01 at f=(1,1) |

f=(0,2) is a frequency the grid resolves more than twenty times over; aliasing
does not explain it.

**Which resolutions are affected is decided entirely by the parity of
`fft_next_fast_size(res x Gmax)`.** At Gmax=11 that is 4, 11 and 12 -- exactly
the three outliers. `res=4` is not special; it is the first resolution that hits
an odd grid at this Gmax.

`(ii + n/2) % n` is bit-identical to the old expression for even `n`, and correct
for odd. Raising the grid to `2*Gmax+1` when the resolution falls short is the
smallest change that brings the index inside; it only ever grows the grid, and
`res >= 3` already clears the condition, so **the default of 8 is untouched.**

| NumBasis | Gmax | grid | buffer | max index | R+T (before) |
|---|---|---|---|---|---|
| 120 | 6 | 12 | 144 | **155** | 1.000064135 |
| 200 | 8 | 16 | 256 | **271** | 0.999991789 |
| 400 | 11 | 24 | 484 | in range | 1.000000000 |

Across the whole probe file those are exactly the two cases where R+T is not 1
to machine precision.

### 6.4 How it was narrowed down

Each step eliminated a candidate rather than confirming a guess. The parity split
survives with the polarization basis off, so it is `DiscretizedEpsilon` itself.
`discretize_cell_probe.c` sums `Pattern_DiscretizeCell`'s area fractions over the
whole grid and finds them exact to 1e-16 for all four shapes, so the
discretisation is innocent. `eps_coeff_probe.cpp` then pulls `Epsilon2` out of
`FMMGetEpsilon_FFT` and compares against the circle's analytic transform -- no
solver, no truncation -- and the **DC coefficient is already wrong**, which puts
it in how the FFT input is filled. Reading the grid traversal finds it.

> **Side finding: `GetEpsilon` does not reflect the discretisation.**
> `Simulation_GetEpsilon` always reconstructs from the analytic
> `Pattern_GetFourierTransform`. Turning on `DiscretizedEpsilon` or
> `SubpixelSmoothing` does not change the return value. Anyone using it to check
> "what structure is S4 actually solving" **gets the wrong answer for the
> discretised path**; the docstring now says so. (The geometric measurements in
> §2 and §3 are on the analytic path and are unaffected.)

### 6.5 After the fix

NumBasis 400, Normal basis:

| res | grid | T before | T after | R+T before | R+T after |
|---|---|---|---|---|---|
| 4 | 45 | 0.871123790 | 0.839006291 | 1.007129241 | 0.999999857 |
| 11 | 125 | 0.843362460 | 0.838541693 | 1.000000000 | 1.000000001 |
| 12 | 135 | 0.840526096 | 0.838529945 | 0.997556894 | 1.000000000 |
| **8 (default)** | 90 | 0.838446275 | 0.838446275 | 1.000000000 | 1.000000000 |

**Every even-grid resolution is bit-identical**, and the odd grids have moved
into the sequence the even ones were tracing. The discretisation-only series
becomes monotone in resolution for the first time:

```
before  0.7956 0.8022 [0.8238] 0.8076 0.8074 ... [0.8128] [0.8125] 0.8093
after   0.7956 0.8022  0.8063  0.8076 0.8074 ...  0.8092   0.8092  0.8093
```

The `res=2` over-read is gone too -- R+T becomes exactly 1 at nb=200.

### 6.6 What this fix does not address

- **`res=2` and `res=3` still alias.** Alias-free coefficients need
  `ngrid > 2·fmax = 4·Gmax`, which resolutions 2 and 3 do not reach. That is the
  consequence of asking for a coarse grid, not a memory error. Overriding the
  resolution the caller asked for is a policy decision rather than a bug fix, so
  it was not done. **All coefficients are resolved only from resolution 4 up.**
- An aliased epsilon matrix is not Hermitian, so R+T can miss 1 by around 1e-6 at
  `res=2` (nb=120 does). That is correct behaviour.

### 6.7 It bites at the default settings too

This is not specific to `res=4`. Which resolution is affected is the parity of
`fft_next_fast_size(res × Gmax)`, and **the default of 8 also comes out odd at
certain Gmax** -- 28, 46, 76, 77, 78.

`Gmax=28` comes from `NumBasis=2453`, an entirely realistic size for a 2D
convergence study. At that setting, with no options touched:

| | DC coefficient error |
|---|---|
| before | **8.869e-03** (the worst coefficient has the wrong sign, at f=(0,−2)) |
| after | **8.882e-16** |

### 6.8 It bites with discretisation off -- the Normal basis uses the same grid

`fmm_PolBasisNV.cpp` also builds a `resolution × Gmax` grid and uses the same
shift. So **`PolarizationDecomposition` alone, without `DiscretizedEpsilon`** --
the path this project's own results run on -- hits the odd grids. NumBasis 400,
analytic epsilon, Normal basis, before the fix:

| res | grid | T |
|---|---|---|
| 3 | 36 | 0.840743440 |
| **4** | **45** | **0.838345383** ← |
| 5 | 60 | 0.840517893 |
| 6 | 72 | 0.840467889 |
| 8 (default) | 90 | 0.840432204 |
| 10 | 120 | 0.840402568 |
| **11** | **125** | **0.839465226** ← |
| **12** | **135** | **0.839532391** ← |
| 13 | 144 | 0.840388073 |

The even ones decrease smoothly; the three odd ones sit 1e-3 to 2e-3 off.

> Incidentally established: **the reference value `0.842114450429` is tied to
> `DiscretizationResolution=8`.** The Normal basis's vector field is built on
> this grid too, so changing the resolution moves the value (at NumBasis 200:
> res 2 → 0.841528, res 8 → 0.842114, res 16 → 0.842174). Reproduction is not at
> risk, since the existing record was taken at the default, but the number is not
> determined by "analytic geometry + Normal" alone.

### 6.9 Direct evidence of the non-determinism

The same calculation gives different answers in different processes. Pre-fix
build, NumBasis 400, Normal, discretised:

| | res=11 T | res=12 T |
|---|---|---|
| run as a standalone script | 0.843362460 | 0.840526096 |
| run in the harness, after other probes | 0.841428300 | 0.841612082 |
| **after the fix (both)** | **0.838541693** | **0.838529945** |

Because the contents of the uninitialised row and column depend on what was
allocated before. After the fix both give the same answer.

### 6.10 Side finding: the Jones basis's R+T excursion was the grid too

The `R+T = 1.0000017` excursion of the `Jones` basis, left unexplained in
MODERNIZE.md §5.2 and in an earlier draft of this report, is **not unexplained
but an under-resolved vector-field grid.** NumBasis 200, analytic epsilon:

| res | R+T − 1 |
|---|---|
| 2 | −2.16e-05 |
| 3 | −3.47e-06 |
| **8 (default)** | **+1.68e-06** |
| 9 | +1.51e-06 |
| 16 | +3.78e-07 |

It converges to 0 as the resolution rises. The Normal basis is already at 1.1e-15
at res=8 under the same conditions (circles and ellipses have a closed-form route
to the normal field and depend on the grid less). **Anyone planning to use Jones
should raise `DiscretizationResolution`.**

Branch `fix/fmm-discretized-grid`, commit `19a0c8c`. Five files, 16 sites (the
shift) plus 5 (the grid lower bound).

---

## 7. The combined build -- this is the new oracle

`modernize/master-fixed` = `modernize/master` + the harness + every fix, merged
without conflict. Compared against the 2016 oracle across all 70 probes:

```
same: 57   differs: 13   failed: 0
```

**failed is 0** -- the three `Clone` probes no longer die. All thirteen
differences are intended:

| probe | why it differs |
|---|---|
| `region_circle_vs_ellipse` | the oracle's axis-aligned ellipse is wrong (`f020606`) |
| `fields_on_grid` | array axis order convention (§2.2 b) |
| `amplitudes_incident` | array length convention (§2.2 c) |
| `output_files` | the transposition above changes `%.14g` by one character |
| `dc_crossing_regions`, `crossing_detection`, `crossing_detection_1d`, `crossing_detection_thin` | overlapping regions are now refused |
| `periodic_image_overlap` | a region meeting a periodic image is now refused |
| `dc_rotated_offcentre_nesting` | the rotation-centre typo is fixed |
| `discretization_grid_parity`, `discretization_grid_overread` | the FMM grid fix |
| `module_entry_points` | `SolveInParallel` and `NewSpectrumSampler` now work |

`two_layer_interface` is in the agreeing column: the oracle was right about a
plane interface all along, and §9 brought this build back to the same numbers.

### 7.1 Installs, and a second build (Python 3.12 / Accelerate)

The environment the experiments run in is Python 3.12, which cannot import a
3.11 extension. It was
rebuilt for 3.12, linked against **macOS Accelerate and S4's bundled kiss_fft**
instead of the conda libraries, so it has no external dependencies at all
(`Accelerate.framework` and `libSystem.B.dylib` only).

| Path | Python | BLAS/FFT | Wiring |
|---|---|---|---|
| a 3.12 install | 3.12 | Accelerate + kiss_fft | one line, a `.pth` in that environment's site-packages |
| a 3.11 install | 3.11 | OpenBLAS + FFTW | `PYTHONPATH` |

**Swapping BLAS and FFT does not change the physics.** Full comparison of the two
builds: `same: 68  differs: 1  failed: 0`. The only difference is 2-3 bytes in
the text file sizes under `output_files` (`%.14g` digit counts).

Running that comparison needs the two builds under different interpreters, hence
`--python-a` / `--python-b` in `compare.py`.

**Added later -- this build procedure was being done by hand.** It was written up
here in prose and never committed as a script, so reproducing it meant reading
prose and retyping four Makefile overrides. It is now `build.sh <dir> accelerate`.

Retyping it also lost one of them. `LIBS` carries `BOOST_LIBS`, which points into the
conda environment, so the rebuilt extension came out with an `@rpath` entry on
`libboost_serialization` -- **a link into the 3.11 conda env, from the
build whose whole purpose is to depend on nothing outside the system.** No source
under `S4/` mentions boost; the flag is vestigial. Cleared along with FFTW3, the
extension is back to `Accelerate.framework` and `libSystem` alone.

The conda env `s4` is untouched: the only thing added is a single `.pth` file on
that side, and the `.so` lives outside the env.

**The reference numbers are unchanged.** On the combined build:

| NumBasis | oracle | combined build | rel |
|---|---|---|---|
| 200 | 0.842114450429 | 0.842114450429 | 6.0e-14 |
| 400 | 0.840432203521 | 0.840432203521 | 5.3e-14 |
| 600 | 0.840907273718 | 0.840907273718 | 1.7e-13 |
| 800 | 0.840475065201 | 0.840475065200 | 6.8e-13 |

The four bases and the off-diagonal complex tensor (`T=0.866227334716`) agree
as well.

---

## 8. Cross-check against an independent implementation

Everything up to here was **S4 against S4**. The 2016 oracle and the 2018 build
share nearly all of their code, so agreement between them means the port did not
change the physics -- **it is not evidence that either one is right.** The only
checks that did not depend on S4's solver were the DC coefficient (geometry
alone) and the Fresnel slab (no pattern at all).

A separate driver fills that gap. It is not in this repository: it reads its
cases and its reference results out of that implementation's tree, so committing it here
would commit a path that only resolves on one machine. It runs the eight cases in
the shared case directory (period 2.5 against wavelength 1.0, so 25 orders propagate and
the comparison is **per order** rather than one number near unity) and diffs S4's
per-order transmission against the other implementation's. The geometry is the non-overlapping
polygons the clipper already resolved, handed identically to both sides.

### 8.1 Two things this cross-check caught

Neither would have been found without running the cross-check it was designed
for.

**(a) A false positive in the crossing test.** S4 refused a case whose two
polygons share **exactly zero** area under a 1200x1200 point-in-polygon test, so
the claim in §5.3 that a refusal always has a witness was wrong. The cause is a
tolerance of zero: clipper output abuts along a shared edge, and the orientation
determinant that should be zero comes back as 4e-18 against terms of size 0.35.
`35a37d7` judges it relative to the terms that formed it, and fixes the two conic
tangency tests for the same reason. Detection is unchanged -- overlaps 11/11,
legitimate 10/10, thin 10/10 -- and `crossing_abutting_clipper_output` pins the
case.

**(b) `LatticeTruncation` had no effect.** The cases use
`NumBasis = (2*fto+1)^2` so that S4's rectangular truncation reproduces
the other implementation's order set exactly, and the counts did not match: 285 instead of 289
at fto 8. `G_select` runs inside `S4_Simulation_New`, before the option can be
given, and nothing re-selects afterwards; the 2016 oracle does it correctly,
which is why it stayed hidden. It is not a small difference -- the rectangular
set at 289 is a 17x17 block where circular truncation discards the corners and
reaches to |m| = 9 -- and **the counts match by chance for some values** (441,
1369), so comparing basis sizes will not reveal it. `a720745` re-runs the
selection when the option changes, passing back the originally requested count,
since `S->n_G` has already been reduced once.

So the first run in §8.2 was rectangular truncation against circular truncation.

### 8.2 Results -- the totals agree, the per-order numbers do not yet

Re-run after the truncation fix, so both sides now use the same rectangular order
set. S4 was pushed to fto 18 (basis 1369).

| Case | other (fto 20) | S4 (fto 18) | total T difference | worst per-order |
|---|---|---|---|---|
| `tilt` (slanted edge) | 0.955030347 | 0.954862527 | **1.7e-04** | 1.00e-02 |
| `overlap` (axis-aligned) | 0.803735331 | 0.802767375 | 9.7e-04 | 2.26e-03 |

**Removing the truncation confound did not change the conclusion** -- the worst
per-order went 1.01e-2 → 1.00e-2 and 2.03e-3 → 2.26e-3, which is to say nowhere.
That defect was real but it is not the cause of this disagreement. (Incidentally
confirmed: after the fix S4's values agree with the oracle's digit for digit --
`overlap` at fto 12 gives 0.802491290 on both -- because the oracle used
rectangular truncation from the start.)

The worst per-order difference **decreases monotonically** with fto: `tilt`
1.99e-2 → 1.00e-2, `overlap` 5.30e-3 → 2.26e-3, with no plateau. That is the
shape of **both being right and neither being converged**, not of one being
wrong.

That `tilt` agrees to 6e-6 in total while individual orders differ by 1e-2 means
**power is being redistributed between orders**, which is exactly what happens
when the Fourier factorisation converges slowly across a slanted boundary.

Which of the three the other implementation paths S4 sits closest to also separates (median over
fto≥10):

| path | median per-order difference from S4 |
|---|---|
| `vector` (block painting) | **7.5e-03** |
| `boundary` | 9.6e-03 |
| `add` (plain superposition) | 5.0e-02 |

Being clearly unlike `add` is independent confirmation that S4 does not simply
sum overlapping contributions.

### 8.3 Narrowing the remaining disagreement

The 9.7e-4 on `overlap` was the sharpest question, since two RCWA
implementations have no business disagreeing on axis-aligned geometry. Three
steps.

**The geometry and the transform are innocent.** `GetEpsilon` FFT'd on a 256x256
grid, against the polygon's analytic transform (divergence-theorem form, checked
to 1e-9 by direct integration): DC agrees to 7.4e-20 and everything inside the
basis to **1.2e-16**.

**The factorisation is not the cause.** The same structure solved four ways --
none, Default, Normal, Jones -- converges on ~0.8028 in all four, within 4e-4 of
each other, and none of them heads for the other implementation's 0.80374.

**The disagreement scales with contrast**: 5.6e-05 at mean eps 1.223824,
1.7e-04 at 1.319464, 9.7e-04 at 1.478032. Which looked like the answer, and
§8.4 explains why it is not.

### 8.4 Correction -- contrast scaling cannot separate the two. 1D can, and they **agree**

The argument in §8.3(3) -- "it scales with contrast, therefore it is a
convergence problem" -- **does not hold as an argument.** An inconsistent
factorisation does not converge slowly; it **converges to the wrong limit** (this
is why Li's rules exist), and that error also scales with contrast. Contrast
scaling cannot distinguish "same limit, different rate" from "different limits".

Distinguishing them means reaching the limit. That is out of reach in 2D and
**within reach in 1D.**

**(1) 1D grating, TE -- agreement at machine precision**

A grating uniform in y (width 0.9 / period 2.5, n=1.6, thickness 0.9, λ=1.0):

| fto | S4 | the other implementation | difference |
|---|---|---|---|
| 4 | 0.9423222188 | 0.9423222188 | 3.7e-15 |
| 14 | 0.9471538644 | 0.9471538644 | 1.3e-13 |

Stack, orders, normalisation, mode solution, S-matrix -- **the machinery is
innocent on both sides.** S4's 1D grating modes, its 2D cell representation and
all three bases give the same value.

**(2) 1D grating, TM -- here they separate**

The same structure differs by 6.8e-4 at fto 14 in p polarisation. A 1D mode has
only `2*fto+1` orders, so very high truncations are cheap. Pushed:

| fto | orders | S4, no decomposition | vs FW | S4 Normal | vs FW |
|---|---|---|---|---|---|
| 14 | 29 | **0.9686005098** | **2.5e-11** | 0.9679199951 | 6.8e-04 |
| 60 | 121 | 0.9686112093 | 1.1e-05 | 0.9684556699 | 1.5e-04 |
| 480 | 961 | 0.9686122569 | 1.2e-05 | 0.9685928319 | 7.7e-06 |

**With `PolarizationDecomposition` off, S4 agrees with the other implementation to 2.5e-11 at
fto 14.** And the `Normal` basis **converges to the same value (0.9686122), but
needs fto 480.**

So **the two implementations are solving the same physics, and their standard
factorisations agree to machine precision.** The disagreement is the convergence
rate of `PolarizationBasis='Normal'` -- a modelling choice, not a defect.

**(3) Back in 2D**

Pushing `overlap` without decomposition to fto 24 (2401 orders) still leaves
0.802942, 7.9e-4 low, rising very slowly (~5e-5 per two fto steps, and not
monotonically). **In 2D the limit is not reachable at any affordable
truncation.** The 1D result explains the same structure's behaviour, so the
remaining 2D difference is likely the same thing -- likely, not proven.

> **Practical conclusion, with its scope corrected**: the observation above is
> **specific to 1D gratings.** In 1D, S4 already has an exact closed-form
> factorisation, and turning on `Normal` **replaces** it with a numerically
> sampled normal field. So this is not a defect but **a misuse: turning Normal on
> in 1D.**
>
> **It must not be generalised to polygon or rectangle geometry.** For curved
> geometry it is the other way round. A cylinder benchmark (period 1, R=0.30,
> n=3.5):
>
> | NumBasis | none | Default | Normal |
> |---|---|---|---|
> | 100 | 0.778451 | 0.834498 | 0.840874 |
> | 400 | 0.809741 | 0.838449 | 0.840432 |
> | 800 | 0.818876 | 0.838923 | 0.840475 |
>
> No decomposition is still 2e-2 away at NumBasis 800 and still climbing. Normal
> is already on its plateau at 100. **For curved geometry, Normal should be on.**
>
> Which is better for 2D polygons (`overlap`) is **not known** -- at fto 18 the
> three methods sit within 4e-4 of each other and none has converged, so there is
> no basis for ranking them.

> **Also worth recording**: the reference series here (Normal, analytic geometry)
> runs 0.842114 / 0.840432 / 0.840907 / 0.840475 over NumBasis 200-800, i.e. it
> **oscillates over a range of 1.7e-3.** It is not a converged value.
> Reproduction is not at risk (same settings, same value) but it should not be
> quoted as "the exact value".

### 8.5 What is established and what is not

**Established**

- Geometry: all eight cases have a mean permittivity exact to between 0 and
  4.4e-16, so S4 holds the structure it was handed
- Total transmission: agreement with an independent implementation to **6e-6** on
  the slanted case
- Direction of convergence: the per-order difference decreases monotonically with
  fto
- No difference from the fixes: oracle median 9.0e-3 vs fixed 9.6e-3 -- **this
  disagreement is not caused by my changes**

**Not established**

- **Per-order agreement is not proven.** It is at the 1e-2 level at fto 18 and
  still falling, but settling it needs fto 24-30 and the cost grows as fto³ (one
  case at fto 18 already takes 135 s).
- On `overlap` (axis-aligned), S4's total is 8.8e-4 below the other implementation's converged
  value and rises by only 4e-5 per two fto steps. Whether it closes is unverified
  -- **the sharpest open question here.**
- The high-contrast cases (`c_19`, `c_22`) cannot be judged, because
  **it has not converged on them either.**

---

## 9. The simplest structure returned nothing

Found by asking a question the harness had never asked: what happens when a
caller does something ordinary that no probe does? Every probe here builds three
layers or more. A plane interface has two.

```python
S.AddLayer(Name='top', Thickness=0, Material='air')
S.AddLayer(Name='bot', Thickness=0, Material='glass')
```

| | T | R | R+T |
|---|---|---|---|
| Fresnel, n = 1 to 1.5 at normal incidence | 0.96 | 0.04 | 1 |
| 2016 oracle | 0.96 | 0.04 | 1 |
| this build, before the fix | **0** | **0** | **0** |

No error. A third layer of any thickness, including zero, gives the right answer,
so the whole defect lives in the two-layer case.

`R+T = 0` is worth pausing on. Energy conservation is the check that catches
almost everything, and it is satisfied here in the sense that nothing came out
because nothing went in -- the incident amplitude itself had been discarded. Only
comparing against an outside answer catches that, which is what §8 was about, and
what the closed form gives for free here.

### 9.1 Where it goes

`SolveAll` moves the two boundary columns to the right-hand side, writing each
into the block row it belongs to:

```c
Copy(n4, t1, 1, &ab[n2], 1);                      /* from a0        */
Copy(n4, t1, 1, &ab[(nlayers-1)*n4-n2], 1);       /* from b_{m-1}   */
```

`(nlayers-1)*n4-n2` equals `n2` when `nlayers` is 2, so the second write lands on
the first. Under front incidence the second contribution is zero, so what it
overwrote the incident amplitude with was zero.

`GetAmplitudes` shows it without any inference: in the top layer the forward
amplitude is the input and the backward amplitude is exactly 0. Nothing scattered
because nothing was there to scatter.

`ConserveMemory=True` selects `SolveInterior` instead, which is unaffected and
returns the Fresnel answer. That the two solve paths disagreed is what localised
the defect to `SolveAll` in one step.

### 9.2 The fix, and what it cannot disturb

Both contributions belong to that block row when it is shared, so they are
summed; for three layers or more the destinations are distinct and the copy is
what is wanted. The guard is on `nlayers == 2` alone.

| | |
|---|---|
| 69 probes on three or more layers, before vs after | `same: 69  differs: 0  failed: 0` |
| two layers vs closed-form Fresnel, 2 indices x 4 angles x 2 polarisations | worst 3.3e-16 |
| both solve paths, same case | agree to the last digit |

### 9.3 The harness gap this exposed

The 62-probe sweep of §2 concluded that the port introduced no physical
difference, and that conclusion stands -- this defect is in `rcwa.cpp`, which the
port did not touch, and it is present in the unmodified port and absent from the
2016 oracle, so it arrived from upstream between the two.

But the sweep could not have found it. Every probe was written by someone who
knew what an S4 script looks like, and an S4 script has an incident half-space, a
structure, and an exit half-space. Two layers looked like a degenerate case not
worth a probe. It is the commonest structure in optics.

`validation/diagnostics.py` now asserts the interface against the closed form
rather than against another build, because on this question neither build is an
authority.

---

## 10. Three things the caller could not find out

Not defects in the physics; gaps in what the binding will tell you about itself.

**Which build.** There is no `__version__`, so an importer cannot tell what it
picked up -- and every build of every commit would have reported the same "1.1"
that upstream declares. It is now the release plus the commit the binary was
built from, with `__release__` and `__build__` separately:

```
>>> S4.__version__
'1.1+e474ef5'          # .dirty appended when S4/ had uncommitted changes
```

An experiment whose result cannot be traced back to the source that produced it
is not reproducible, and for a fork that changes solver behaviour that is not a
nicety.

**What the options are.** `SetOptions` was write-only. Nothing confirmed that an
option had been understood, which is exactly how `LatticeTruncation` being
ignored (§8.1b) stayed hidden -- the only symptom was an order count, and
circular truncation returns the requested count for some values. `GetOptions`
returns the keys `SetOptions` accepts, so its result feeds straight back in, and
reports `NumBasis` beside `NumBasisRequested`:

```
LatticeTruncation='Circular'         NumBasis 285  (requested 289)
LatticeTruncation='Parallelogramic'  NumBasis 289  (requested 289)
```

**What the methods are called.** Every `SetRegion*` docstring named its
pre-rename Lua name -- `SetRegionCircle` documented itself as
`SetLayerPatternCircle` -- so `help()` pointed at an API that has not existed for
years. All 31 are rewritten against the real keywords, and the two sharp edges
say so where they will actually be read: `GetEpsilon` does not reflect
`DiscretizedEpsilon`, and `GetFieldsOnGrid` indexes `[x][y]` where the 2016
release indexed `[y][x]`.

Branch `fix/introspection`, commit `07b0ef2`.

---

## 11. The caller's side -- accepting both keyword spellings

`cb47b74` prefixed the basic types with `S4_` and ran that rename over the Python
keyword-argument strings as well, so `Layer=` / `Material=` became `S4_Layer=` /
`S4_Material=`. Nothing else about the API changed.

The cost is borne entirely by the caller.

- A script written before it raises `TypeError` on a build after it, and the
  reverse.
- The message names **the keyword the caller wrote**, not the one now wanted, so
  the reader has to already know the rename happened.
- The two builds are **indistinguishable from Python.** `validation/s4compat.py`
  tells them apart by making a call and watching it fail.

Both spellings now work. A bare keyword is renamed in **a copy of the keyword
dict**, so nothing the caller owns is touched, and the prefixed spelling wins if
both are somehow given.

The parse goes through `S4_ParseKeywords` instead of
`PyArg_ParseTupleAndKeywords` so that the temporary dict's release lives in one
place rather than on every return path of thirty-three methods.

```python
S.SetRegionCircle(Layer='pat',    Material='rod', Center=(0,0), Radius=0.2)
S.SetRegionCircle(S4_Layer='pat', S4_Material='rod', Center=(0,0), Radius=0.2)
# T = 0.968659190421806 -- the two lines give the same number
```

`diagnostics.py` asserts it, **against the raw module rather than through the
shim** -- the shim exists to hide this difference and would hide the test with it
-- and asserts that the two spellings give the same number, not merely that
neither raises.

No numeric regression (`same: 69  differs: 0`). Branch
`fix/keyword-compatibility`, commit `5b9345c`.

---

## 12. Baselines, and the order the fixes have to go in

Local branch names are not recorded here: they move, and twice in one afternoon
they made this section wrong. Two things about the layout are durable and worth
keeping.

**The baselines**, because the isolation experiments in §2.2 and §9 build
against them and cannot be reproduced without them:

| | commit | what it is |
|---|---|---|
| `baseline/upstream` | `7fd00a2` | victorliu/S4 master, untouched |
| `baseline/port-unfixed` | `21302a5` | the port, before any of these fixes |

A pull request to upstream should be cut from `baseline/upstream`. Everything
described in this report sits on top of it.

**The order**, because three of the fixes depend on each other for reasons in
the code rather than in the history:

- The crossing test must sit above the rotation-centre fix. Its interior-point
  search calls `shape_contains_point`, so with the typo present it refuses
  **legitimate** nested pairs whose parent is rotated and off the diagonal --
  measured, not assumed (§5.7).
- The periodic-image check, and then the message detail, each attach indices and
  pair information to the check below them, so they follow it.
- The keyword shim touches every parse path in `main_python.c`, so it goes last
  of the binding changes.

Everything else is independent and merges in any order, and each fix carries
only C changes: the tests live in the harness, so the same probes run before and
after from any point.

---

## 13. What could not be established

- **Aliasing at `res=2` and `res=3` remains** (§6.6). It is the consequence of
  asking for a coarse grid rather than a memory error, and overriding the
  resolution the caller gave is a policy decision, so it was left alone.
- **False negatives in the crossing test -- reduced, not zero.** §5.6 changed the
  failure condition. It used to miss overlaps thin relative to **the shapes**;
  now it can miss overlaps small relative to **the bounding-box intersection**.
  Such shapes usually have properly crossing boundaries and are caught by the
  segment test, so the real gap is much narrower, but it is not zero. Only the
  absence of false positives is guaranteed.
- **The meaning of `GetSMatrixDeterminant`.** It is the determinant of a 4n×4n
  matrix and underflows at any practical basis size. It agrees with the oracle at
  NumBasis 1/5/9, but at 120 the two builds produce different approximations of
  zero. Do not rely on it for anything beyond sign tracking in a mode search.
- **The other frontends were not built.** The Lua message table was updated, but
  Lua, MPI, Matlab and R were not compiled.
- **The CHOLMOD path is unverified** (disabled in this build).
- **Built with `-mcpu=apple-m1`.** Not portable to another machine as is.
