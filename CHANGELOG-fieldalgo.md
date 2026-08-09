# s4-fieldalgo

A fork of [victorliu/S4](https://github.com/victorliu/S4), whose last release is
1.1 and whose last commit is June 2018. This line carries fixes to defects in
that code, and versions itself separately: reporting 1.1 would be a claim about
behaviour that is no longer true.

The import name is unchanged. `import S4` works, and so does every script
written against either keyword spelling (`Layer=` and `S4_Layer=`).

```python
>>> import S4
>>> S4.__version__          # this line's release, plus the commit it was built from
'0.1.0+902b081'
>>> S4.__upstream__
'victorliu/S4 1.1 (7fd00a2)'
```

`upstream/ChangeLog` is upstream's and is left alone; this file covers only what
this line changed.

---

## 0.1.0

Fourteen defects fixed, all of them present in upstream master. `UPSTREAM.md`
gives each one its smallest reproduction; `validation/REPORT.md` gives what was
measured. `validation/diagnostics.py` asserts the lot in about a minute.

### Structures that used to be accepted and are now refused

These are the changes most likely to affect an existing script. In every case
the structure could not be represented and was previously solved as something
else, silently.

- Two regions in a layer whose boundaries cross. Nesting and disjointness are
  unaffected; a region tangent to another, or sharing a boundary point, still
  passes.
- A region overlapping a periodic image of itself or of another region. A region
  may still cross the unit cell boundary; only overlap once tiled is refused.
- Overlapping intervals in a 1D layer, including nesting, which 1D has never
  been able to represent.
- A repeated layer name. `SetLayer` remains the way to redefine a layer.
- A lattice basis spanning no cell -- parallel vectors, or a 1D period of zero.
  A zero second vector (how a 1D lattice is written) and a negatively oriented
  basis are both still accepted.

The refusals name the layer and the region, counting from 1 in the order the
regions were added.

### Results that change

- **A structure of exactly two layers returned nothing at all.** A plane
  interface gave T = 0, R = 0, R+T = 0. It now matches Fresnel to 3.3e-16.
  Three layers or more were never affected.
- The discretised-epsilon grid dropped a row and a column into uninitialised
  memory whenever the grid came out odd, which the default resolution reaches at
  some basis sizes, and read past its buffer at resolution 2. Even-grid results
  are bit-identical; odd-grid ones move into the sequence the even ones traced.
- `shape_contains_point` rotated about the wrong centre, so a region nested
  inside a rotated, off-diagonal-centre parent was recorded as sitting on the
  background.
- The ellipse tangent code dropped a term, indexed by the loop count rather than
  the loop variable, and read past a four-element array. It reached every
  polarization basis, and produced answers that varied between runs.
- `LatticeTruncation` had no effect from Python: the basis was chosen before the
  option could be set.
- `Clone()` produced an unusable clone and damaged the simulation it came from.

### New

- `GetOptions()` returns the options in force, with the keys `SetOptions`
  accepts, plus `NumBasis` beside `NumBasisRequested` -- the difference between
  them is how a wrong truncation shows itself.
- `__version__`, `__release__`, `__build__`, `__upstream__`. The build id is the
  commit, with `.dirty` when the working tree had uncommitted changes under
  `S4/`.
- `SolveInParallel` does something. `NewSpectrumSampler` no longer raises
  `SystemError`.
- `SetVerbosity` can be called at all; its parse check was inverted.

### Fixed with no change to any result

- `GetEpsilon` returned a value with an exception set, so the real error was
  replaced by `SystemError`.
- The rectangle normal used `-.1` where `-1.` was meant; normalisation absorbed
  it.
- Diagnostic messages were truncated mid-character for long or non-ASCII layer
  names, and lost the part that said what was wrong.
- Every `SetRegion*` docstring named its pre-rename Lua function. All 31 now
  describe the real keywords, and note the two sharp edges: `GetEpsilon` does
  not reflect `DiscretizedEpsilon`, and `GetFieldsOnGrid` indexes `[x][y]` where
  the 2016 release indexed `[y][x]`.

### Known limits

`validation/REPORT.md` §13 has the full list. The ones worth knowing before
using this for anything: resolutions 2 and 3 still alias; the crossing test has
no false positives but its false-negative rate is not zero; per-order agreement
with an independent implementation is not established beyond fto 18; and the
reference convergence series is tied to `DiscretizationResolution=8`.

---

## Releasing

The version lives in `version.txt` and nowhere else -- `gensetup.py.sh` reads it
into both the C source and the distribution metadata, so the two cannot
disagree.

1. Add the entry here.
2. Edit `version.txt`.
3. Commit, then `git tag -a v<version>`.
4. Rebuild the installs: `bash validation/build.sh <dir> [conda|accelerate]`.
   The build stamps itself with the commit, so tag first or the stamp will name
   the commit before the tag.
5. `cd validation && python diagnostics.py`.

Installs are kept rather than overwritten. Each build lands in its own slot
under `s4-builds/`, named `<version>+<commit>-py<X.Y>-<flavour>`, and the path
given on the command line becomes a symlink to it:

```
s4-builds/
    0.1.0+b8b24fb-py3.11-conda/       S4...so, BUILD-INFO.txt
    0.1.0+b8b24fb-py3.12-accelerate/
    INDEX.txt                         one line per slot
s4-fixed-py312 -> s4-builds/0.1.0+b8b24fb-py3.12-accelerate
```

Switching back to an earlier build is repointing the symlink; nothing that
consumes the path has to change. They are about 400 KB each, and each is the
only thing that can answer what produced a number from an experiment already
run -- rebuilding from the tag gives the same source, but Makefile.local, the
compiler and the environment's libraries all sit outside the tag.

The archive sits beside the checkout rather than inside it, and `builds` in the
repository root is a symlink to it. `git clean -xdf` is what one reaches for
when a build is broken, which is the worst possible moment to lose every
earlier build; this way it removes the link and leaves the archive.

Set `S4_BUILD_ARCHIVE` to keep them somewhere else.

## Publishing

The public branch is one commit on top of upstream's `7fd00a2`, carrying the
whole tree; the development history stays local. Nothing about it is secret --
it simply is not what a reader of this fork needs, and squashing keeps the
public log to the thing that is: what changed against upstream, and how to
check it.

```sh
NEW=$(git commit-tree 'dev^{tree}' -p 7fd00a2 -m "$MESSAGE")
git branch -f master $NEW            # dev keeps the full history
git tag -f -a v<version> $NEW
git push origin master v<version>
```

`git diff dev master` must be empty afterwards: the squash changes the history,
never the tree.

Push tags by name. `--tags` would also push `baseline/port-unfixed`, which is
not an ancestor of the public commit and would drag its own history along with
it -- that one is wanted (the isolation experiments build against it, and it
carries the ellipse fix the reports cite), but as a deliberate choice rather
than a side effect.

## Version numbers

Bump the **minor** when solver output changes for a structure that used to be
accepted, or when something that used to be accepted is refused -- both mean an
existing script can get a different answer or stop running. Bump the **patch**
for anything that cannot change a result: messages, documentation, tests, and
fixes to paths that were already unreachable.

It is 0.x because the surface is still moving. `GetEpsilon` not reflecting the
discretised path is unresolved, and there is no one-call geometry check yet.
