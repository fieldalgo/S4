#!/bin/sh
# Build libS4 under AddressSanitizer and run validation/discretize_probe.c on it.
#
#   sh validation/asan.sh [NumBasis] [DiscretizationResolution]
#
# Defaults to 120 2, the smallest configuration where the discretised-epsilon
# grid comes out at exactly 2*Gmax and the coefficient read runs off the end of
# the FFT buffer.
#
# This leaves build/ holding sanitized objects; run validation/build.sh
# afterwards to get an ordinary library back.
set -e

NB=${1:-120}
RES=${2:-2}
ROOT=$(cd "$(dirname "$0")/.." && pwd)
cd "$ROOT"

CONDA=${CONDA_PREFIX:-/opt/miniconda3/envs/s4}
export PATH=$CONDA/bin:$PATH

test -f Makefile.local || cp Makefile.local.example Makefile.local

SAN="-fsanitize=address -fno-omit-frame-pointer -g -O1"

rm -rf build
make -f Makefile.local objdir
make -f Makefile.local -j4 \
	CFLAGS="-Wall -m64 -mcpu=apple-m1 -mtune=native -fPIC $SAN" \
	OPTFLAGS="-O1 $SAN"

cc -I S4 $SAN -o /tmp/discretize_probe validation/discretize_probe.c \
	build/libS4.a -lm -lstdc++ \
	-L$CONDA/lib -lopenblas -lfftw3 -Wl,-rpath,$CONDA/lib

echo "--- running NumBasis=$NB resolution=$RES ---"
ASAN_OPTIONS=detect_leaks=0 /tmp/discretize_probe "$NB" "$RES"
