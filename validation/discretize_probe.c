/* Drive the discretised-epsilon path from C, so it can be run under a sanitizer.
 *
 * The Fourier coefficients of a discretised layer are read out of the FFT
 * output (fmm_FFT.cpp, fmm_kottke.cpp, fmm_PolBasis*.cpp all do this) with
 *
 *     f[i] = G[2*a+i] - G[2*b+i];  if(f[i] < 0){ f[i] += ngrid[i]; }
 *     ... Fto[f[1] + f[0]*ngrid[1]]
 *
 * f ranges over +/- 2*Gmax, and the wrap only repairs the negative half.  The
 * read is therefore in bounds only while ngrid > 2*Gmax.  The grid is
 * fft_next_fast_size(resolution * Gmax), which at resolution 2 is exactly
 * 2*Gmax whenever that product already factors into 2s, 3s and 5s -- and 2 is
 * the documented minimum resolution.  Then f[0] can equal ngrid[0] and the read
 * runs past the end of a buffer of ngrid[0]*ngrid[1] elements:
 *
 *     NumBasis 120 -> Gmax 6, resolution 2 -> ngrid 12, ng2 144, index up to 155
 *     NumBasis 200 -> Gmax 8, resolution 2 -> ngrid 16, ng2 256, index up to 271
 *
 * Both of those configurations break energy conservation (R+T = 1.000064135 and
 * 0.999991789 against 1.0 to machine precision everywhere else), which is what
 * this exists to attribute.
 *
 * Build and run: see validation/asan.sh.
 *
 * STATUS: the sanitized binary does not start on this machine.  It deadlocks
 * before main() -- no output at all, including a write to stderr on the first
 * line -- whether it is linked against the conda OpenBLAS or against
 * Accelerate.  The conda environment ships its own libclang_rt.asan, and the
 * first attempt died on a version-mismatch symbol from it, so that was the
 * obvious suspect; removing it from the link did not help, and the cause is
 * still unknown.  The over-read above is therefore established by the index
 * arithmetic and by the energy-conservation break, not by a sanitizer report.
 * This file is committed so that whoever picks the sanitizer up next starts
 * from a driver that already sets the case up.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "S4.h"
#include "S4_internal.h"

int main(int argc, char *argv[]) {
	const int nb  = (argc > 1) ? atoi(argv[1]) : 120;
	const int res = (argc > 2) ? atoi(argv[2]) : 2;

	const S4_real Lr[4] = {1, 0, 0, 1};
	S4_Simulation *S = S4_Simulation_New(Lr, (unsigned)nb, NULL);

	S4_real eps_air[2] = {1.0, 0.0};
	S4_real eps_rod[2] = {3.5 * 3.5, 0.0};
	S4_MaterialID air = S4_Simulation_SetMaterial(
		S, -1, "air", S4_MATERIAL_TYPE_SCALAR_COMPLEX, eps_air);
	S4_MaterialID rod = S4_Simulation_SetMaterial(
		S, -1, "rod", S4_MATERIAL_TYPE_SCALAR_COMPLEX, eps_rod);

	S4_real t0 = 0.0, t1 = 0.15;
	S4_LayerID top = S4_Simulation_SetLayer(S, -1, "top", &t0, -1, air);
	S4_LayerID pat = S4_Simulation_SetLayer(S, -1, "pat", &t1, -1, air);
	S4_LayerID bot = S4_Simulation_SetLayer(S, -1, "bot", &t0, -1, air);

	{
		const double centre[2] = {0.0, 0.0};
		Simulation_AddLayerPatternCircle(S, &S->layer[pat], (int)rod, centre, 0.30);
	}

	S->options.use_discretized_epsilon = 1;
	S->options.resolution = res;
	S->options.use_polarization_basis = 1;
	S->options.use_normal_vector_basis = 1;

	{
		double angle[2] = {0.0, 0.0};
		double pol_s[2] = {1.0, 0.0};   /* magnitude, phase */
		double pol_p[2] = {0.0, 0.0};
		Simulation_MakeExcitationPlanewave(S, angle, pol_s, pol_p, 0);
	}
	S->omega[0] = 2.0 * M_PI / 1.6;
	S->omega[1] = 0.0;

	{
		double power_top[4], power_bot[4];
		double off = 0.0;
		int r0 = S4_Simulation_GetPowerFlux(S, top, &off, power_top);
		int r1 = S4_Simulation_GetPowerFlux(S, bot, &off, power_bot);
		if (0 != r0 || 0 != r1) {
			fprintf(stderr, "solve failed: %d %d\n", r0, r1);
			return 1;
		}
		printf("nb=%d res=%d  T=%.9f  R+T=%.9f\n", nb, res,
		       power_bot[0] / power_top[0],
		       (power_bot[0] - power_top[1]) / power_top[0]);
	}

	S4_Simulation_Destroy(S);
	return 0;
}
