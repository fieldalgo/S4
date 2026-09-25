# Defects found in victorliu/S4 at master (7fd00a2)

Written to be sent upstream. Everything here reproduces on unmodified master with
the Python extension built on macOS arm64 / clang 17 / Python 3.11; nothing
depends on this fork's changes except where a fix is offered.

Each item gives the smallest input that shows it, what actually happens, and why.
The reproductions are Python because that is the frontend under test, but ten of
the fourteen are in the C core and reach the Lua frontend identically.

Severity is about what a user gets, not how hard it was to find:

- **silent** — a wrong number comes back with no error, warning or message
- **crash** — the process dies
- **loud** — an error is raised, but the wrong one, or where none was warranted

| # | Defect | File | Severity |
|---|---|---|---|
| 1 | `S4_Simulation_Clone` leaves the clone unusable and damages the source | `S4/S4.cpp` | crash |
| 2 | `shape_contains_point` rotates about the wrong centre | `S4/pattern/pattern.c` | silent |
| 3 | Regions whose boundaries cross are accepted and mis-solved | `S4/pattern/pattern.c` | silent |
| 4 | A region overlapping its own periodic image is accepted | `S4/pattern/pattern.c` | silent |
| 5 | The FMM grid shift is not a permutation for an odd grid | `S4/fmm/*.cpp` | silent |
| 6 | The FMM coefficient read can index past its buffer | `S4/fmm/*.cpp` | silent |
| 7 | `LatticeTruncation` has no effect from Python | `S4/main_python.c` | silent |
| 8 | `SolveInParallel` does nothing at all | `S4/main_python.c` | silent |
| 9 | `NewSpectrumSampler`'s format string is malformed | `S4/main_python.c` | loud |
| 10 | A structure of exactly two layers returns nothing at all | `S4/rcwa.cpp` | silent |
| 11 | `SetVerbosity` cannot be called at all | `S4/main_python.c` | loud |
| 12 | A repeated layer name is accepted | `S4/S4.cpp` | silent |
| 13 | A lattice basis with no cell is accepted | `S4/S4.cpp` | silent |
| 14 | The ellipse tangent code drops a term, indexes by the loop count, and reads out of bounds | `S4/pattern/pattern.c` | silent |
| 15 | `S4_Simulation_SetMessageHandler` has no `return` statement | `S4/S4.cpp` | crash |

Two more that are real but change no result, listed at the end.

Number 10 is the one to read first if you read only one. It is the simplest
structure the code can express, and it comes back empty.

---

## 1. `S4_Simulation_Clone` — the clone is unusable and the source is damaged

`S4.cpp`, `S4_Simulation_Clone`.

```python
S = S4.New(Lattice=((1,0),(0,1)), NumBasis=20)
S.SetMaterial(Name='air', Epsilon=1.0)
S.AddLayer(Name='t', Thickness=0, S4_Material='air')
S.AddLayer(Name='b', Thickness=0, S4_Material='air')
S.SetExcitationPlanewave(IncidenceAngles=(0,0), sAmplitude=1, pAmplitude=0)
S.SetFrequency(1/1.6)
S.GetPowerFlux(S4_Layer='t')
C = S.Clone()
C.GetPowerFlux(S4_Layer='t')      # SIGSEGV
```

`memcpy(T, S, sizeof *S)` copies the counts along with everything else, and the
materials and layers are then added through `S4_Simulation_SetMaterial`/`SetLayer`
with an id of `-1`, which appends at the current count. The entries land at
`[n, 2n)` and `[0, n)` keeps whatever `malloc` returned; every by-name lookup walks
from zero and `strcmp`s against an uninitialised pointer.

Three more faults sit behind it:

- The setter is handed `M->type`, the internal encoding (0 scalar, 1 tensor), where
  it expects an `S4_MATERIAL_TYPE_*` constant — those start at 2, so the value
  always reaches the switch's `default`, which discards the material just appended
  and copies no epsilon. A clone surviving the first fault would still have had no
  materials.
- `S4_Simulation_SetLayer` calls `Simulation_DestroySolution` on its target. With
  `T->solution` still aliasing `S->solution`, adding the clone's first layer frees
  the **source's** solution and leaves `S` pointing at it.
- `G`, `kx`, the vector-field dump prefix and polygon vertex arrays stay aliased for
  the clone's whole life, so destroying both simulations frees each twice.

The 2016 revision wrote into slots `[i]` directly and did not have this.

Fix offered: `6bbac7f` on this fork.

---

## 2. `shape_contains_point` rotates about the wrong centre

`pattern.c`:

```c
x[0] = (x_[0] - s->center[0]) * ca + (x_[1] - s->center[1]) * sa;
x[1] = (x_[0] - s->center[1]) *-sa + (x_[1] - s->center[1]) * ca;
/*              ^^^^^^^^^ center[0] */
```

The query point is displaced by `sin(angle) * (center[1] - center[0])`, a constant
independent of the query. It vanishes for an unrotated region and on the diagonal
`cx == cy`, which is why tests that rotate about the origin or offset without
rotating never see it. It also needs two or more regions, because with one the
containment tree does no comparisons.

The predicate feeds `pattern_get_containment_tree`, so a region genuinely nested
inside a rotated, off-diagonally centred parent is recorded as sitting on the
background and its Fourier contribution is taken against the wrong permittivity.
That is a different structure, not a loss of accuracy.

Measured on the DC Fourier coefficient, which is exact whatever the truncation —
an eps-4 circle of radius 0.04 at the centre of an eps-9 rectangle of half-widths
(0.40, 0.08) rotated 45° about (0.15, −0.15), background eps 1, cell 1×1:

| | mean permittivity |
|---|---|
| exact | 1.998867258770 |
| S4 | 2.039079644740 |

2.039079644740 is exactly the value for the circle sitting on the background. With
`Angle=0` and everything else identical, S4 is exact to 4e-16.

Fix offered: `cd9f4a9`.

---

## 3. Regions whose boundaries cross are accepted and mis-solved

The containment tree records one parent per region and the transform subtracts each
region from exactly that parent, so it can express nesting and disjointness and
nothing else. A pair whose boundaries cross is neither, and the tree silently picks
one of the two readings.

Two eps-4 and eps-9 squares of half-width 0.17 centred at x = ∓0.08, background
eps 1, cell 1×1:

| | mean permittivity |
|---|---|
| exact (later region wins) | 2.0880 |
| S4 | 1.9248 |

The eps-4 square's exposed area contributes nothing.

`pattern_get_containment_tree` means to reject this. Two things stop it. The pair
loop advances and tests `i` in both positions:

```c
for(j = i+1; i < nshapes; ++i){
    if(shapes_intersect(&shapes[i], &shapes[j])){ ... }
}
```

so the outer loop is finished after one pass and the inner one compares shapes
against a fixed partner, itself included. And `shapes_intersect` is a stub — every
branch of its type switch falls through to `return 0`, with the one real test
commented out.

1D is worse. A 1D region is a rectangle with a zero second half-width, so its area
is zero and `shape_contains_point` rejects every point; the tree leaves all of them
parented to the background and the transform adds each interval's full
contribution. Overlapping intervals double-count the shared part whether the
overlap is partial or total, so 1D has no working notion of nesting either — a
region written inside another reads 8.3 where 5.9 is exact.

Fix offered: `35a37d7`, which implements `shapes_intersect` (exact segment crossing
for polygons, the quadratic in the conic's normalised frame for polygon-conic, a
sign-change sweep for conic-conic, plus an interior-overlap test for pairs whose
boundaries only meet along collinear edges), repairs the loop, and adds the 1D
interval check. Note that the orientation tests need a tolerance relative to the
products that form them: clipper output abuts along shared edges and the
determinant that should be zero comes back as a few ulps.

The 1D refusal wants an error code of its own rather than the 2D one, whose remedy
is "nested or disjoint" and is false there: nesting is not available in 1D at all.
The diagnostic also has to name the region, which costs nothing — the containment
tree already encodes the index in its return value, and both checks run before the
sort, so it is the order the caller added the regions in.

---

## 4. A region overlapping its own periodic image is accepted

A region may leave the unit cell — the transform is taken over the shape as given,
which is right for a periodic structure exactly while the shape and its repeats are
disjoint. Once they overlap the shared area is counted once per copy.

A circle of radius 0.7 in a 1×1 cell, eps 9 on eps 1:

```
mean permittivity reported: 13.3
maximum possible (cell entirely rod): 9
```

The value is outside the range the materials allow, and no error is raised.

The pair need not be one region twice. A small region near the cell edge landing
inside a large region's image double-counts identically, and a 1D lattice makes
both forms unreachable for any area-based test: a 1D region has zero height, so it
intersects nothing. Two intervals at ∓0.45 with half-width 0.1 meet across the cell
boundary, and a single interval of half-width 0.6 runs into its own repeats; both
are accepted.

Fix offered: `bb8c277`, extended in `3d4a036` to report both indices of the pair
and to reduce the 1D separation modulo the period. Touching an image is not an
overlap, which has to be kept: a rectangle spanning the cell exactly is how a 1D
grating is written inside a 2D lattice, and an interval exactly as wide as the
period tiles it correctly.

---

## 5. The FMM grid shift is not a permutation for an odd grid

`fmm_FFT.cpp:90` and fifteen more sites across `fmm_kottke.cpp`,
`fmm_PolBasisNV.cpp`, `fmm_PolBasisJones.cpp`, `fmm_PolBasisVL.cpp`:

```c
si = ii >= n/2 ? ii - n/2 : ii + n/2;
```

With `n` even this is a rotation by `n/2`. With `n` odd it is not a bijection —
`n/2` truncates, so for `n = 5` the indices map to `{2,3,0,1,2}`: slot 2 is written
twice and slot `n-1` is never written. In two dimensions a whole row and a whole
column of the sample grid keep whatever `S4_malloc` returned, and every Fourier
coefficient is computed from them.

The damage is not confined to high frequencies. At NumBasis 400 with
`DiscretizationResolution` 4 (grid 45), the DC coefficient itself is off by 4.4e-2
where an even grid gets it to 1.3e-14, and the worst coefficient error sits at
f = (0,2) — a frequency the grid resolves twenty times over — reading −0.229 where
the exact value is +0.042.

Which resolutions are affected depends only on whether
`fft_next_fast_size(resolution * Gmax)` is odd. **The default resolution of 8 is
not safe**: Gmax = 28, 46, 76, 77 and 78 all produce odd grids, and Gmax = 28 is
NumBasis 2453. At that setting, with no option touched, the DC coefficient is off
by 8.9e-3.

Because the affected slots are uninitialised, the answer can change between runs.
On the 2016 revision the same input gave three different results in twenty runs
(twelve `T=0`, five `T=0.8639`, three `NaN`), with LAPACK reporting
`ZGEEV parameter number 4 had an illegal value`.

`(ii + n/2) % n` is the same map for even `n`, bit for bit, and a proper rotation
for odd `n`.

Fix offered: `19a0c8c`.

---

## 6. The FMM coefficient read can index past its buffer

Same files:

```c
f = G_a - G_b;  if(f < 0){ f += ngrid; }  ...  Fto[f[1] + f[0]*ngrid[1]]
```

`f` spans ±2·Gmax and only the negative half is wrapped, so the index is valid only
while `ngrid > 2*Gmax`. The grid is `fft_next_fast_size(resolution * Gmax)`, which
at resolution 2 — the documented minimum — is exactly `2*Gmax` whenever that product
already factors into 2s, 3s and 5s:

| NumBasis | Gmax | grid | buffer | index reaches | R+T |
|---|---|---|---|---|---|
| 120 | 6 | 12 | 144 | 155 | 1.000064135 |
| 200 | 8 | 16 | 256 | 271 | 0.999991789 |

Those two were the only configurations in a 65-case sweep whose R+T was not 1 to
machine precision.

Fix offered: raising the grid to `2*Gmax+1` when the multiplier does not reach it,
in `19a0c8c`. It can only increase the grid and only at resolution 2.

---

## 7. `LatticeTruncation` has no effect from Python

```python
S = S4.New(Lattice=((2.5,0),(0,2.5)), NumBasis=289)
S.SetOptions(LatticeTruncation='Parallelogramic')
# ... solve ...
len(S.GetBasisSet())     # 285, and the set is circular
```

`G_select` runs inside `S4_Simulation_New` from `S->options.lattice_truncation`,
which is still the default there because `New` takes only a lattice and a count.
Nothing re-selects afterwards, and `S4Sim_SetLatticeTruncation` is commented out of
the method table, so from Python there is no way to obtain a rectangular basis at
all. The 2016 revision honours the option.

The difference is not cosmetic. At NumBasis 289 the rectangular set is the 17×17
block of orders; circular truncation drops the block's corners and reaches out to
|m| = 9. The **count** sometimes matches by coincidence — 441 and 1369 both come out
exact — while the set does not, so comparing basis sizes shows nothing wrong.

Fix offered: `a720745`, re-running the selection when the option changes. It has to
be handed the count the caller originally asked for, because `S->n_G` has already
been reduced by the first selection.

---

## 8. `SolveInParallel` does nothing

`main_python.c`:

```c
static PyObject *S4_SolveInParallel(PyObject *Self, PyObject *args, PyObject *kwds)
{
	static char *kwlist[] = {"S4_Layer", "Simulations", NULL};
	const char *layerName;
	//S4_solve_in
	Py_RETURN_NONE;
}
```

It never calls the parser. `S4.SolveInParallel()` with no arguments at all returns
`None`. The caller is told the layers have been solved when nothing has happened.
Results stay correct only because the layer is solved again on demand.

Fix offered: `1a9047d`, which parses the arguments, checks each entry is a
simulation, resolves the layer by name and solves it. The work is serial — the Lua
frontend threads it through helpers bound to a `lua_State` that this binding has no
equivalent of — and the docstring says so.

---

## 9. `NewSpectrumSampler`'s format string is malformed

```c
PyArg_ParseTupleAndKeywords(args, kwds, "dd|i|d|d|d|O!:SpectrumSampler_New", ...)
```

`|` appears six times; CPython accepts it once. Any call supplying an optional
argument raises `SystemError: Invalid format string (| specified twice)`, so the
sampler works only when given its two required arguments and nothing else.

`"dd|idddO!"` fixes it. Fix offered: `1a9047d`.

---

## 10. A structure of exactly two layers returns nothing at all

A plane interface -- two layers, nothing between them -- is the simplest structure
S4 can express and the only one with a closed-form answer. It comes back empty.

```python
S = S4.New(Lattice=((1,0),(0,1)), NumBasis=9)
S.SetMaterial(Name='air',   Epsilon=1.0)
S.SetMaterial(Name='glass', Epsilon=2.25)          # n = 1.5
S.AddLayer(Name='top', Thickness=0, S4_Material='air')
S.AddLayer(Name='bot', Thickness=0, S4_Material='glass')
S.SetExcitationPlanewave(IncidenceAngles=(0,0), sAmplitude=1, pAmplitude=0)
S.SetFrequency(1/1.6)
```

| | T | R | R+T |
|---|---|---|---|
| Fresnel | 0.96 | 0.04 | 1 |
| S4 | **0** | **0** | **0** |

No error, no warning. Adding a third layer of any thickness, even zero, gives the
right answer, so the defect is confined to two.

`SolveAll` moves the two boundary columns to the right-hand side by writing each
into the block row it belongs to (`rcwa.cpp`, "Prepare the RHS"):

```c
Copy(n4, t1, 1, &ab[n2], 1);                      /* from a0 */
Copy(n4, t1, 1, &ab[(nlayers-1)*n4-n2], 1);       /* from b_{m-1} */
```

`(nlayers-1)*n4-n2` is `n2` when `nlayers` is 2, so the second write lands on the
first. Under front incidence the second contribution is zero, so the incident
amplitude is overwritten with zero and the structure is solved as though nothing
had been sent into it. `GetAmplitudes` shows it directly: the forward amplitude
in the top layer is the input, and the backward amplitude is exactly 0.

`SolveInterior` -- the path `ConserveMemory` selects -- is not affected and gives
the Fresnel answer, which is the quickest confirmation:

```python
S.SetOptions(ConserveMemory=True)     # T = 0.96, R = 0.04
```

Fix offered: sum the two contributions when they share a block row, which they do
only when `nlayers == 2`. Every one of 69 existing probes on three or more layers
is bit-identical afterwards, and the two-layer case then matches Fresnel to
3.3e-16 across two indices, four angles and both polarisations.

---

## 11. `SetVerbosity` cannot be called at all

```c
if(PyArg_ParseTupleAndKeywords(args, kwds, "i:SetVerbosity", kwlist, &level))
    return NULL;
```

The test is inverted. A call that parses correctly returns NULL without setting
an exception, and CPython reports that as

```
SystemError: <method 'SetVerbosity' of 'S4.Simulation' objects> returned NULL
without setting an exception
```

A call that does not parse falls through to use an uninitialised `level`, on top
of the exception `PyArg_ParseTupleAndKeywords` has already set -- and writes it
into `options.verbosity` whenever the stack happens to hold a value in [0, 9].

So the option is unreachable from Python. Every other parse site in the file is
written with the negation; this is the only one that is not.

Fix offered: add the `!`.

---

## 12. A repeated layer name is accepted

```python
S.AddLayer(Name='slab', Thickness=0.2, S4_Material='glass')
S.AddLayer(Name='slab', Thickness=0.5, S4_Material='glass')
```

Both layers enter the stack and both change the answer -- transmission goes from
0.529 to 0.397 -- while every lookup by name resolves to the first, so the second
cannot be addressed, patterned or measured. `GetPowerFlux(Layer='slab')` reports
the first.

`Simulation_InitSolution` reserves error code 12 for exactly this and the check
that would return it is commented out. It is also written inside the non-copy
branch, so two layer copies sharing a name would pass even with it enabled.

Redefining a layer is what `SetLayer` is for, so a repeated name from `AddLayer`
is a mistake rather than an idiom.

Fix offered: enable the check for every layer, and name both indices in the
message.

---

## 13. A lattice basis with no cell is accepted

```python
S4.New(Lattice=((1,0),(2,0)), NumBasis=9)     # parallel vectors
S4.New(Lattice=0.0, NumBasis=9)               # 1D period of zero
```

Both span no unit cell, so the reciprocal lattice does not exist. What comes back
is `NaN` for every result, and on stderr -- not as an exception --

```
** On entry to  ZGEEV parameter number  4 had an illegal value
Layer eigensystem returned info = -4
```

`S4_Lattice_Reciprocate` already detects this and returns 1 for parallel vectors
and 2 for both zero. `S4_Simulation_New` discards the return value, and the code
that would have reported it is commented out immediately below the call.

A zero second vector is not degenerate -- that is how a 1D lattice is written --
and neither is a negatively oriented basis; both must keep working.

Fix offered: check in `Simulation_InitSolution`, the one gate both frontends pass
through, so the error is raised where it can be reported.

---

## 14. The ellipse branch of the tangent-cross code has three faults

`shape_get_tangent_cross_segment` and its `_tri` twin build the constraints the
normal vector field is solved from, so they run for every polarization basis --
`Normal`, `Jones` and the default -- whenever a layer holds an ellipse. The
ellipse branch of both (`pattern.c:588` and `:670`) is wrong in three ways at
once.

```c
u1 = d0[0]*-sa;
u[0] = d0[0]*ca; u[1] = u1 * ratio;          /* (1) */

c = intersection_circle_segment(s->vtab.ellipse.halfwidth[0], p1, u, isect, NULL);
for(i = 0; i < c; ++i){
    double L;
    isect[2*c+1] *= iratio;                  /* (2) */
    L = hypot(isect[2*c+0],isect[2*c+1]);    /* (2) */
    if(0 == L){ L = 1; }
    *cross += (u[0]*isect[2*c+0] + u1*isect[2*c+1]) / L;   /* (2), (3) */
}
```

**(1)** The segment direction is rotated into the ellipse frame using only its
x component. `u[0]` is missing `+ d0[1]*sa` and `u[1]` is missing `d0[1]*ca`, so
the rotation is effectively discarded. The position `p0p` two lines above is
rotated correctly, which is what the direction should have matched.

**(2)** Every subscript inside the loop uses `c`, the number of intersections,
where `i`, the loop variable, is meant. `isect` is `double[4]` -- two points --
so with two intersections `isect[2*c+0]` is `isect[4]`, **one element past the
end of the array**. The loop also revisits the same slot on every pass instead
of walking the intersections.

**(3)** The dot product uses `u1`, which is the y component before the `ratio`
scaling, rather than `u[1]`.

The effect is visible in an identity that needs no reference implementation: an
ellipse with equal half-widths **is** a circle, and rotating a circle changes
nothing.

```python
S.SetOptions(PolarizationDecomposition=True, PolarizationBasis='Normal')
# ... same layer, same off-centre position, once as each shape
```

| shape | T |
|---|---|
| circle, R = 0.30 | 0.807627216212 |
| ellipse (0.30, 0.30) at 0 deg | 0.782190513141 |
| ellipse (0.30, 0.30) at 25, 45, 70, 90 deg | 0.765815666420 (all four) |

Two things are wrong at once. Writing the circle as an ellipse moves the answer
by 4.2e-2, and the ellipse is not even rotation-invariant, though a circle must
be: 0 degrees sits 1.6e-2 away from every other angle. Under `Jones` and the
default basis the 0-degree case comes back `NaN`.

Zero is the special one because `sa` vanishes there, so the term missing from
`u[0]` costs nothing and only `u[1]` is wrong; at any other angle both
components are. Why the remaining angles agree with each other has not been
traced.

Rotation itself does reach the calculation -- an ellipse of unequal half-widths
varies by 0.15 in T across these angles, before the fix as well as after -- so
this is not the rotation being discarded. The companion identity says the same:
an ellipse rotated 90 degrees against the same ellipse with its half-widths
swapped agrees to 8e-15 even unfixed. The geometry is fine; the direction
handed to the intersection routine is not.

### Scope

Isolated by reverting `f020606` alone on top of an otherwise identical build, so
nothing else moves. Of 70 probes, 2 change.

- **Only with the polarization decomposition on.** With it off the two builds are
  bit-identical, including for ellipses -- the tangent code is reached only
  through `Pattern_GenerateFlowField`, which only the `PolBasis*` formulations
  call.
- **Only ellipses.** Circles, rectangles and polygons are bit-identical.
- **Not only the circle case.** An ellipse of unequal half-widths is wrong too:
  (0.34, 0.18) moves by 1.5e-2 at 0 degrees and 4.0e-2 at 90.
- **Non-deterministic at 0 degrees.** Fault (2) reads `isect[4]`, one element
  past a four-element array, so what comes back depends on the stack. The same
  command run eight times in eight processes, ellipse (0.34, 0.18) at 0 degrees:

  | run | `Normal` | `Jones` | `Default` |
  |---|---|---|---|
  | 1, 2, 4, 8 | 0.9190417985 | `nan` | `nan` |
  | 3, 6 | 0.9190417985 | 0.4917745978 | `nan` |
  | 5 | 0.9190417985 | 0.8435506600 | `nan` |
  | 7 | 0.9190417985 | 0.4917745978 | 1 |

  with `ZGEEV parameter number 4 had an illegal value` on stderr each time. The
  fixed build gives the same three numbers in all eight runs.

- **The basis stops mattering.** Away from 0 degrees the three polarization
  bases agree with each other to 3e-8, where after the fix they differ by up to
  1.8e-2 as different factorisations should. The decomposition is being computed
  from a field that no longer distinguishes them.

### What the fix is verified against

Afterwards the circle and the equal-half-width ellipse agree to 1.0e-13 across
`Normal`, `Jones` and the default basis at seven angles. That identity is exact
and needs no reference implementation, and the circle takes an entirely separate
branch of the same function, so it is an independent route to the same answer.

A second route agrees for that case: approximating the ellipse by an N-gon --
which goes through `intersection_polygon_segment` instead -- converges onto the
fixed value (5.5e-3 at 24 sides down to 3.9e-5 at 384) and stays 4.2e-2 away
from the unfixed one at every N.

For an ellipse of **unequal** half-widths there is no exact identity, and the
polygon route cannot stand in for one. Refining the polygon does not converge to
the ellipse under the normal-vector basis; it plateaus at 1.7e-3 and 6.2e-3 and
stops moving between 384 and 768 sides. With the decomposition off the same
refinement converges cleanly, gaining a factor of four per doubling
(4.0e-4 to 1.6e-6), so the geometry and its transform are right and the residue
is the normal field itself: a faceted polygon's normals keep their corners
however fine it gets, while an ellipse's come from the smooth analytic tangent.

So the fix is confirmed exactly where an identity exists, and elsewhere is only
shown to be far closer to an independent construction than the code it replaces.

---

## 15. `S4_Simulation_SetMessageHandler` has no `return` statement

`S4.cpp`, `S4_Simulation_SetMessageHandler`.

It is declared to return an `S4_message_handler` and falls off the end:

```c
S4_message_handler S4_Simulation_SetMessageHandler(
	S4_Simulation *S, S4_message_handler handler, void *data
){
	if(NULL == S){ return NULL; }
	S->msg = handler;
	S->msgdata = data;
}
```

Flowing off the end of a non-void function is undefined behaviour in C++, and the
compiler is entitled to assume the path is never taken. Unlike 1-14 this is latent
on upstream master: neither frontend calls the function, so it is reachable only
through the C API. This fork's Python binding installs a handler in `S4Sim_new`
(`main_python.c`), which puts it on the path of every `New`:

```python
S = S4.New(Lattice=((1,0),(0,1)), NumBasis=9)   # *** stack smashing detected ***: terminated
```

Measured on Ubuntu 24.04 x86_64, gcc 13.3, Python 3.12, OpenBLAS 0.3.26, with the
library built from `Makefile.local.example` at `-O3`: the process aborts before
`New` returns. Rebuilt at `-O1 -g`, gdb shows the call from `S4Sim_new` arriving in
`S4_Simulation_GetLattice` -- the next function in the object file -- with the
handler's address as `Lr`, and its `memcpy` segfaults. The build says so up front:
`S4.cpp:397: warning: control reaches end of non-void function [-Wreturn-type]`.
The macOS arm64 / clang builds this fork was validated on cannot have hit it --
every diagnostic calls `New` -- which undefined behaviour permits, and is why it
went unnoticed.

Fix offered: return the handler being replaced, which is what the return type
implies (the `signal()` convention). With it, `validation/diagnostics.py` passes
25/25 on the Linux build above.

---

## Two that change no result

**`GetEpsilon` returns a value with an exception set.** It calls
`HandleSolutionErrorCode` on a non-zero return and then falls through to build a
complex from the untouched output buffer. CPython replaces the real error with
`SystemError: ... returned a result with an exception set`, so the caller sees
neither the message nor a usable number. Every other entry point in the file
already returns `NULL` there. Fix in `35a37d7`.

**`shape_get_normal` writes `-.1` for `-1.`**, twice, in the RECTANGLE branch. It
changes no answer: the function normalises before returning, so the vector comes
back unit-length pointing the same way and the only trace is one ulp on two of
twelve components. Worth correcting anyway, because the only consumer today
happens to normalise again and a caller trusting the documented unit length would
be handed a vector ten times too short on two of the four faces. Fix in `470ef1b`.

---

## How these were found

A sweep that calls every binding entry point and compares the numbers, digit for
digit, against the 2016 revision — which shares almost all of this code, so
agreement between them means the two agree, not that either is right. Geometry is
checked on the DC Fourier coefficient, which is exact whatever the truncation and
so separates "the wrong structure was built" from "the solve has not converged".

Items 1, 7, 8 and 9 are visible only by calling the entry point. Items 2, 3 and 4
need the DC coefficient. Items 5 and 6 needed reading the FMM after the DC
coefficient showed the epsilon matrix was wrong at a frequency the grid resolves
easily.

The harness is `validation/` in this fork; `validation/REPORT.md` carries the
measurements and the reasoning for each fix.
