/* Does Pattern_DiscretizeCell tile the unit cell?
 *
 * The discretised-epsilon path asks Pattern_DiscretizeCell, for every pixel of
 * an nu x nv grid, what fraction of that pixel each region covers.  Summed over
 * the whole grid those fractions have to reproduce each region's area exactly,
 * whatever nu and nv are -- that is the one property the FFT downstream relies
 * on, and it is what fixes the DC Fourier coefficient.
 *
 * No FFT and no solver here: just the tiling.
 *
 * Usage: discretize_cell_probe [shape]
 *   shape is one of circle (default), rectangle, ellipse, polygon.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pattern/pattern.h"

static const double RAD = 0.30;

int main(int argc, char *argv[]) {
	const char *what = (argc > 1) ? argv[1] : "circle";
	const double L[4] = {1, 0, 0, 1};
	Pattern p;
	shape s;
	int parent[1];
	double exact;
	double polyv[8] = {-0.22, -0.14, 0.24, -0.16, 0.19, 0.21, -0.14, 0.23};
	int nu;

	memset(&s, 0, sizeof(s));
	s.center[0] = 0.0;
	s.center[1] = 0.0;
	s.angle = 0.0;
	if (0 == strcmp(what, "rectangle")) {
		s.type = RECTANGLE;
		s.vtab.rectangle.halfwidth[0] = 0.31;
		s.vtab.rectangle.halfwidth[1] = 0.17;
		exact = 4 * 0.31 * 0.17;
	} else if (0 == strcmp(what, "ellipse")) {
		s.type = ELLIPSE;
		s.vtab.ellipse.halfwidth[0] = 0.33;
		s.vtab.ellipse.halfwidth[1] = 0.19;
		exact = M_PI * 0.33 * 0.19;
	} else if (0 == strcmp(what, "polygon")) {
		s.type = POLYGON;
		s.vtab.polygon.n_vertices = 4;
		s.vtab.polygon.vertex = polyv;
		exact = 0.0;
		{ /* shoelace */
			int i, j;
			for (j = 0, i = 3; j < 4; i = j++) {
				exact += polyv[2 * i + 0] * polyv[2 * j + 1]
				       - polyv[2 * j + 0] * polyv[2 * i + 1];
			}
			exact = fabs(exact) / 2;
		}
	} else {
		s.type = CIRCLE;
		s.vtab.circle.radius = RAD;
		exact = M_PI * RAD * RAD;
	}

	parent[0] = -1;
	p.nshapes = 1;
	p.shapes = &s;
	p.parent = parent;

	printf("%s: exact area = %.12f\n", what, exact);
	printf("  %4s %4s %16s %14s\n", "nu", "par", "sum of fractions", "error");
	for (nu = 6; nu <= 65; ++nu) {
		double total = 0.0;
		int iu, iv;
		for (iu = 0; iu < nu; ++iu) {
			for (iv = 0; iv < nu; ++iv) {
				double value[2] = {0, 0};
				Pattern_DiscretizeCell(&p, L, nu, nu, iu, iv, value);
				total += value[1];
			}
		}
		total /= (double)(nu * nu);
		if (nu <= 20 || fabs(total - exact) > 1e-9) {
			printf("  %4d %4s %16.12f %14.3e%s\n", nu, (nu % 2) ? "odd" : "even",
			       total, fabs(total - exact),
			       (fabs(total - exact) > 1e-9) ? "   <<<" : "");
		}
	}
	return 0;
}
