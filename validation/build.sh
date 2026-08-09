#!/bin/sh
# Build libS4 and install the Python extension into an isolated directory.
#
#   sh validation/build.sh /tmp/s4-somewhere              # 3.11, conda OpenBLAS + FFTW
#   sh validation/build.sh ~/project/fieldalgo/s4-fixed-py312 accelerate
#
# Never run `make S4_pyext`: its last line is `pip3 install --upgrade ./`, which
# installs into whatever environment happens to be active.  The install below
# goes to --target and touches nothing else.  The guard around it fails the
# build if the shared conda env's extension changes for any reason.
#
# The Makefile does not track header dependencies, so editing a header leaves
# stale objects behind; this removes build/ every time rather than trying to be
# clever about when that matters.
set -e

DEST=${1:?usage: build.sh <install-dir> [conda|accelerate]}
FLAVOR=${2:-conda}
ROOT=$(cd "$(dirname "$0")/.." && pwd)
cd "$ROOT"

# Which interpreter to build for.  S4_PYTHON overrides; otherwise each flavour
# has a default that is used only if it exists, and anything else falls back to
# whatever python3 is on PATH.  Nothing here should require a particular
# machine's layout.
pick_python() {
	if [ -n "$S4_PYTHON" ]; then
		echo "$(cd "$(dirname "$S4_PYTHON")" && pwd)"
	elif [ -x "$1/python" ]; then
		echo "$1"
	else
		command -v python3 >/dev/null || {
			echo "FATAL: no python3 on PATH and $1 does not exist" >&2; exit 1; }
		dirname "$(command -v python3)"
	fi
}

case "$FLAVOR" in
conda)
	# The 3.11 env that also holds the 2016 oracle, so both sides of a
	# comparison see the same OpenBLAS and the same FFTW.
	PYBIN=$(pick_python /opt/miniconda3/envs/s4/bin)
	set --
	;;
accelerate)
	# Accelerate for BLAS/LAPACK and the kiss_fft S4 already ships, so the
	# extension links against system libraries only (Accelerate.framework and
	# libSystem) and can be dropped into an environment that has neither
	# OpenBLAS nor FFTW.  Empty FFTW3_LIB is what selects kiss_fft: the
	# Makefile's `ifdef` is false for an empty value.
	#
	# BOOST_LIBS goes too.  No source under S4/ mentions boost -- the flag is
	# vestigial -- but it is in LIBS, so it ends up as an @rpath entry pointing
	# into the 3.11 conda env, which is exactly the dependency this flavour
	# exists to avoid.
	# No local default here: this flavour exists to depend on nothing outside
	# the system, so it takes whatever python3 is on PATH unless S4_PYTHON says
	# otherwise.  Point S4_PYTHON at the interpreter the install is for --
	# building against the wrong minor version produces an extension that will
	# not import.
	PYBIN=$(pick_python "")
	set -- "BLAS_LIB=-framework Accelerate" "LAPACK_LIB=-framework Accelerate" \
		"FFTW3_INC=" "FFTW3_LIB=" "BOOST_INC=" "BOOST_LIBS="
	;;
*)
	echo "unknown flavor '$FLAVOR' (want conda or accelerate)" >&2
	exit 1
	;;
esac

export PATH=$PYBIN:$PATH

# An extension installed in a shared environment that this build must not
# touch: `make S4_pyext` would pip-install over it.  Set S4_GUARD to point at
# yours, or leave it unset and the check is skipped.
GUARD=${S4_GUARD:-/opt/miniconda3/envs/s4/lib/python3.11/site-packages/S4.cpython-311-darwin.so}
if [ -f "$GUARD" ]; then
	BEFORE=$(shasum -a 256 "$GUARD" | cut -d' ' -f1)
else
	GUARD=""
fi

if [ ! -f Makefile.local ]; then
	cp Makefile.local.example Makefile.local
	echo "Makefile.local created from the template; edit it if this build needs"
	echo "a BLAS, an FFT or a compiler other than the defaults."
fi

rm -rf build
make -f Makefile.local "$@" objdir
make -f Makefile.local "$@" -j4

# Stamp the binary with the source it came from.  An experiment whose result
# cannot be traced to a build is not reproducible, and the working tree is
# checked separately from HEAD because a binary built over uncommitted edits is
# exactly the one that later cannot be explained.
BUILD_ID=$(git rev-parse --short HEAD 2>/dev/null || echo unknown)
if ! git diff --quiet HEAD -- S4 gensetup.py.sh 2>/dev/null; then
	BUILD_ID="$BUILD_ID.dirty"
fi
export S4_BUILD_ID=$BUILD_ID

printf 'printlibs:\n\t@echo $(LIBS)\n' > /tmp/s4-printlibs.mk
LIBSVAL=$(make -f Makefile.local -f /tmp/s4-printlibs.mk "$@" -s printlibs)
sh gensetup.py.sh ./build ./build/libS4.a "$LIBSVAL"

# Installs are kept, not overwritten.  Each one is ~400 KB and each is the only
# thing that can answer "what produced this number" for an experiment already
# run -- rebuilding from the tag gives the same source but not necessarily the
# same binary, since Makefile.local, the compiler and the environment's
# libraries all sit outside the tag.  The path the caller asked for becomes a
# symlink to the slot, so everything that consumes it keeps working.
PYTAG=$("$PYBIN/python" -c 'import sys;print("%d.%d"%sys.version_info[:2])')
FORKVER=$(cat version.txt 2>/dev/null || echo 0.0.0)
ARCHIVE=${S4_BUILD_ARCHIVE:-$(cd "$(dirname "$DEST")" && pwd)/s4-builds}
SLOT="$ARCHIVE/$FORKVER+$BUILD_ID-py$PYTAG-$FLAVOR"

mkdir -p "$ARCHIVE"
rm -rf "$SLOT"
"$PYBIN/pip" install . --no-build-isolation --target "$SLOT" >/dev/null

# Record how this install was made.  The version stamp names the source; it
# cannot name the configuration, and Makefile.local -- which decides the
# compiler flags and the whole link line -- is deliberately outside version
# control because it carries absolute paths for one machine.  Without this, two
# installs stamped with the same commit can differ and nothing says why.
{
	echo "s4-fieldalgo $(cat version.txt 2>/dev/null || echo unknown)+$BUILD_ID"
	echo "flavour     $FLAVOR"
	echo "commit      $(git rev-parse HEAD 2>/dev/null || echo unknown)"
	echo "branch      $(git rev-parse --abbrev-ref HEAD 2>/dev/null || echo unknown)"
	echo "describe    $(git describe --tags --always --dirty 2>/dev/null || echo unknown)"
	echo "built       $(date -u '+%Y-%m-%dT%H:%M:%SZ')"
	echo "host        $(uname -smr)"
	echo "compiler    $(g++ --version 2>/dev/null | head -1)"
	echo "python      $("$PYBIN/python" -c 'import sys;print(sys.version.split()[0])' 2>/dev/null)"
	echo "makefile    sha256 $(shasum -a 256 Makefile.local | cut -d' ' -f1)"
	echo "libs        $LIBSVAL"
} > "$SLOT/BUILD-INFO.txt"

rm -rf "$DEST"
ln -sfn "$SLOT" "$DEST"

# One line per slot, so the archive can be read without opening each manifest.
{
	for d in "$ARCHIVE"/*/; do
		[ -f "$d/BUILD-INFO.txt" ] || continue
		printf '%-44s %s\n' "$(basename "$d")" \
			"$(sed -n 's/^built  *//p' "$d/BUILD-INFO.txt")"
	done
} | sort > "$ARCHIVE/INDEX.txt"

if [ -n "$GUARD" ]; then
	AFTER=$(shasum -a 256 "$GUARD" | cut -d' ' -f1)
	if [ "$BEFORE" != "$AFTER" ]; then
		echo "FATAL: the guarded install changed during this build" >&2
		echo "  before $BEFORE" >&2
		echo "  after  $AFTER" >&2
		exit 1
	fi
fi
echo "built $BUILD_ID -> $SLOT"
echo "  $DEST -> $(basename "$SLOT")"
[ -n "$GUARD" ] && echo "guard intact: $AFTER"
