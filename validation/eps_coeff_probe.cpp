/* Compare the discretised-epsilon Fourier coefficients against the exact ones.
 *
 * FMMGetEpsilon_FFT builds the 2n x 2n Epsilon2 matrix by sampling the layer on
 * a grid and FFTing it.  Entry (i,j) of its upper-left block is the coefficient
 * at G_i - G_j, so the whole difference set can be read out and checked against
 * the analytic transform of a circle,
 *
 *     eps_hat(q) = (eps_rod - eps_bg) * R * J1(2*pi*R*q) / q ,  q = |G|,
 *     eps_hat(0) =  eps_bg + (eps_rod - eps_bg) * pi * R^2 ,
 *
 * with no solve involved and no truncation error on either side.  That isolates
 * which coefficients the grid gets wrong, which is what the run-to-run and
 * grid-parity anomalies come down to.
 *
 * Usage: eps_coeff_probe <NumBasis> <resolution> [--dump]
 *   --dump lists every distinct |G| with its computed and exact coefficient.
 */

#include <complex>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>

#include "S4.h"
#include "S4_internal.h"
#include "fmm/fmm.h"

static const double RAD = 0.30;
static const double EPS_BG = 1.0;
static const double EPS_ROD = 3.5 * 3.5;

/* J1 via the standard library; the coefficient is real for a centred circle. */
static double exact_coeff(double q) {
	if (q <= 0) {
		return EPS_BG + (EPS_ROD - EPS_BG) * M_PI * RAD * RAD;
	}
	return (EPS_ROD - EPS_BG) * RAD * j1(2.0 * M_PI * RAD * q) / q;
}

static int next_fast_size(int n) {
	for (;;) {
		int m = n;
		while (m % 2 == 0) m /= 2;
		while (m % 3 == 0) m /= 3;
		while (m % 5 == 0) m /= 5;
		if (m <= 1) return n;
		++n;
	}
}

int main(int argc, char *argv[]) {
	const int nb  = (argc > 1) ? atoi(argv[1]) : 400;
	const int res = (argc > 2) ? atoi(argv[2]) : 4;
	const int dump = (argc > 3);

	const S4_real Lr[4] = {1, 0, 0, 1};
	S4_Simulation *S = S4_Simulation_New(Lr, (unsigned)nb, NULL);

	S4_real eps_air[2] = {EPS_BG, 0.0};
	S4_real eps_rod[2] = {EPS_ROD, 0.0};
	S4_MaterialID air = S4_Simulation_SetMaterial(
		S, -1, "air", S4_MATERIAL_TYPE_SCALAR_COMPLEX, eps_air);
	S4_MaterialID rod = S4_Simulation_SetMaterial(
		S, -1, "rod", S4_MATERIAL_TYPE_SCALAR_COMPLEX, eps_rod);

	S4_real t1 = 0.15;
	S4_LayerID pat = S4_Simulation_SetLayer(S, -1, "pat", &t1, -1, air);
	{
		const double centre[2] = {0.0, 0.0};
		Simulation_AddLayerPatternCircle(S, &S->layer[pat], (int)rod, centre, RAD);
	}
	S->layer[pat].pattern.parent =
		(int*)malloc(sizeof(int) * S->layer[pat].pattern.nshapes);
	Pattern_GetContainmentTree(&S->layer[pat].pattern);

	S->options.use_discretized_epsilon = 1;
	S->options.resolution = res;
	S->omega[0] = 2.0 * M_PI / 1.6;
	S->omega[1] = 0.0;

	const int n = S->n_G;
	const int n2 = 2 * n;
	std::complex<double> *Epsilon2   = new std::complex<double>[n2 * n2];
	std::complex<double> *Epsilon_inv = new std::complex<double>[n * n];

	FMMGetEpsilon_FFT(S, &S->layer[pat], n, Epsilon2, Epsilon_inv);

	int gmax = 0;
	for (int g = 0; g < n; ++g) {
		if (abs(S->G[2 * g + 0]) > gmax) gmax = abs(S->G[2 * g + 0]);
		if (abs(S->G[2 * g + 1]) > gmax) gmax = abs(S->G[2 * g + 1]);
	}
	/* Mirror the library's grid rule.  Two variants, because the point of this
	 * probe is to compare a build that has the lower bound against one that does
	 * not: report both and let the caller see which the numbers belong to. */
	const int ngrid_nobound = next_fast_size(gmax * res);
	int ngrid_bounded = gmax * res;
	if (ngrid_bounded <= 2 * gmax) ngrid_bounded = 2 * gmax + 1;
	ngrid_bounded = next_fast_size(ngrid_bounded);
	const int ngrid = ngrid_bounded;

	/* Walk every (i,j) pair, key by the difference vector, keep the worst error
	 * seen for each shell of |f|. */
	std::map<int, double> worst;   /* f[0]*f[0]+f[1]*f[1] -> max abs error */
	std::map<int, double> mag;     /* the same shell -> |exact| */
	double worst_all = 0.0;
	int worst_f[2] = {0, 0};
	double worst_got = 0.0, worst_want = 0.0;

	for (int j = 0; j < n; ++j) {
		for (int i = 0; i < n; ++i) {
			const int f0 = S->G[2 * i + 0] - S->G[2 * j + 0];
			const int f1 = S->G[2 * i + 1] - S->G[2 * j + 1];
			const double q = sqrt((double)(f0 * f0 + f1 * f1));
			const double want = exact_coeff(q);
			const std::complex<double> got = Epsilon2[i + j * n2];
			const double err = std::abs(got - std::complex<double>(want, 0.0));
			const int shell = f0 * f0 + f1 * f1;
			if (worst.find(shell) == worst.end() || err > worst[shell]) {
				worst[shell] = err;
				mag[shell] = fabs(want);
			}
			if (err > worst_all) {
				worst_all = err;
				worst_f[0] = f0; worst_f[1] = f1;
				worst_got = got.real(); worst_want = want;
			}
		}
	}

	printf("nb=%d res=%d  Gmax=%d  fmax=%d  grid: unbounded-rule %d (%s), "
	       "bounded-rule %d (%s)\n",
	       nb, res, gmax, 2 * gmax,
	       ngrid_nobound, (ngrid_nobound % 2) ? "ODD" : "even",
	       ngrid_bounded, (ngrid_bounded % 2) ? "ODD" : "even");
	printf("  worst coefficient error %.6e at f=(%d,%d): got %.9f  exact %.9f\n",
	       worst_all, worst_f[0], worst_f[1], worst_got, worst_want);
	printf("  unbounded rule would index %s; bounded rule indexes %s\n",
	       (ngrid_nobound > 2 * gmax) ? "in bounds" : "OUT OF BOUNDS",
	       (ngrid_bounded > 2 * gmax) ? "in bounds" : "OUT OF BOUNDS");
	printf("  note: coefficients are free of aliasing only while ngrid > 2*fmax"
	       " = %d\n", 4 * gmax);

	if (dump) {
		printf("  %8s %14s %14s\n", "|f|^2", "max abs err", "|exact|");
		for (std::map<int, double>::const_iterator it = worst.begin();
		     it != worst.end(); ++it) {
			printf("  %8d %14.6e %14.6e\n", it->first, it->second, mag[it->first]);
		}
	}

	delete [] Epsilon2;
	delete [] Epsilon_inv;
	S4_Simulation_Destroy(S);
	return 0;
}
