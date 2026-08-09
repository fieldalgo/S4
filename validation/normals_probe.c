/* Print the region normals the library computes, to full precision.
 *
 * shape_get_normal is not reachable from the Python binding: the only path that
 * consumes it is the Kottke subpixel-averaging branch of the FMM, which folds
 * the vector into a rotation and never reports it.  A defect in that function
 * is therefore invisible to every Python-level test, so this links against
 * libS4.a and asks it directly.
 *
 * The RECTANGLE branch writes
 *     double sgn = (rx > 0. ? 1. : -.1);
 * where -1. is plainly meant.  Whether that matters at all is what this
 * measures: shape_get_normal normalises before returning, so a vector built
 * from -0.1 comes back with the same direction and unit length as one built
 * from -1., and only the rounding differs.  Run this against a build with and
 * without the correction and diff the output.
 *
 * Build (from the repository root, after `make -f Makefile.local`):
 *     cc -I S4 -o /tmp/normals_probe validation/normals_probe.c build/libS4.a -lm -lstdc++
 */

#include <math.h>
#include <stdio.h>
#include <string.h>

#include "pattern/pattern.h"

static void report(const char *label, const shape *s, double x, double y) {
	double p[2], n[2];
	p[0] = x;
	p[1] = y;
	shape_get_normal(s, p, n);
	printf("%-28s (% .3f,% .3f) -> (% .17g, % .17g)  |n| = %.17g\n",
	       label, x, y, n[0], n[1], hypot(n[0], n[1]));
}

int main(void) {
	shape s;
	double ang;
	int k;

	/* A rectangle at a general angle, so that both branches of the sign test
	 * are reached on all four sides and neither ca nor sa is 0 or 1. */
	for (k = 0; k < 3; ++k) {
		const double angles[3] = {0.0, 30.0, 57.0};
		char label[64];
		ang = angles[k] * M_PI / 180.0;
		memset(&s, 0, sizeof(s));
		s.type = RECTANGLE;
		s.center[0] = 0.05;
		s.center[1] = -0.03;
		s.angle = ang;
		s.vtab.rectangle.halfwidth[0] = 0.30;
		s.vtab.rectangle.halfwidth[1] = 0.12;

		printf("rectangle, angle = %g deg, halfwidths (0.30, 0.12)\n", angles[k]);
		/* one point just outside each face, in the rectangle's own frame */
		{
			const double ca = cos(ang), sa = sin(ang);
			const double faces[4][2] = {
				{ 0.40,  0.00},   /* +x face */
				{-0.40,  0.00},   /* -x face: the branch with the sign typo */
				{ 0.00,  0.20},   /* +y face */
				{ 0.00, -0.20},   /* -y face: the other branch with the typo */
			};
			const char *names[4] = {"+x face", "-x face (sign)",
			                        "+y face", "-y face (sign)"};
			int f;
			for (f = 0; f < 4; ++f) {
				const double gx = s.center[0] + ca * faces[f][0] - sa * faces[f][1];
				const double gy = s.center[1] + sa * faces[f][0] + ca * faces[f][1];
				snprintf(label, sizeof(label), "  %s", names[f]);
				report(label, &s, gx, gy);
			}
		}
		printf("\n");
	}

	/* A circle and an ellipse for contrast: neither branch uses a sign
	 * constant, so these pin that nothing else moved. */
	memset(&s, 0, sizeof(s));
	s.type = CIRCLE;
	s.center[0] = 0.05;
	s.center[1] = -0.03;
	s.vtab.circle.radius = 0.25;
	printf("circle, radius 0.25\n");
	report("  right", &s, 0.35, -0.03);
	report("  upper left", &s, -0.15, 0.14);
	printf("\n");

	memset(&s, 0, sizeof(s));
	s.type = ELLIPSE;
	s.center[0] = 0.05;
	s.center[1] = -0.03;
	s.angle = 30.0 * M_PI / 180.0;
	s.vtab.ellipse.halfwidth[0] = 0.30;
	s.vtab.ellipse.halfwidth[1] = 0.12;
	printf("ellipse, angle 30 deg, halfwidths (0.30, 0.12)\n");
	report("  right", &s, 0.35, -0.03);
	report("  upper left", &s, -0.15, 0.14);

	return 0;
}
