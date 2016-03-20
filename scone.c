#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include "planets.h"
#include "vectors.h"

double ea_from_ma(double ma, double ecc)
{
	double left = 0, right = M_PI * 2.0, mid, val;

	ma = fmod(ma, 2 * M_PI);
	while (right - left > 1e-9) {
		mid = (left + right) / 2.0;
		val = mid - ecc * sin(mid);
		if (val < ma)
			left = mid;
		else
			right = mid;
	}
	return mid;
}

static inline double radius_at_ta(struct planet p, double true_anomaly)
{
	return p.sma * (1 - p.ecc * p.ecc) / (1.0 + p.ecc * cos(true_anomaly));
}

struct vector position_at_ut(struct planet p, double ut)
{
	double n = sqrt(sungp / pow(p.sma, 3));
	double period = 2.0 * M_PI / n;
	double mean_anomaly = 2.0 * M_PI * (ut - epoch) / period;
	double ecc_anomaly = ea_from_ma(mean_anomaly, p.ecc);
	double true_anomaly = 2.0 * atan2(sqrt(1.0 + p.ecc) * sin(ecc_anomaly / 2.0),
					  sqrt(1.0 - p.ecc) * cos(ecc_anomaly / 2.0));
	struct vector v = {radius_at_ta(p, true_anomaly), 0, 0};

	v = rotateXY(v, true_anomaly + p.ape);
	v = rotateYZ(v, p.inc);
	v = rotateXY(v, p.lan);
	return v;
}

struct vector velocity_at_ut(struct planet p, double ut)
{
	double n = sqrt(sungp / pow(p.sma, 3));
	double period = 2.0 * M_PI / n;
	double mean_anomaly = 2.0 * M_PI * (ut - epoch) / period;
	double ecc_anomaly = ea_from_ma(mean_anomaly, p.ecc);
	double true_anomaly = 2.0 * atan2(sqrt(1.0 + p.ecc) * sin(ecc_anomaly / 2.0),
					  sqrt(1.0 - p.ecc) * cos(ecc_anomaly / 2.0));
	/* ṙ = (na/sqrt(1-e²)){-sin ν, e+cos ν, 0} */
	double vscale = n * p.sma / sqrt(1 - p.ecc*p.ecc);
	struct vector v = {-vscale*sin(true_anomaly), vscale*(p.ecc + cos(true_anomaly)), 0};

	v = rotateXY(v, p.ape);
	v = rotateYZ(v, p.inc);
	v = rotateXY(v, p.lan);
	return v;
}

double ve(struct planet p, double r)
{
	return sqrt(2.0 * p.gp / (p.r + r));
}

/* vhe is hyperbolic excess velocity */
double v_from_vhe(struct planet p, double vhe, double r)
{
	/* v² = ve² + vhe² */
	return sqrt(2.0 * p.gp / (p.r + r) + vhe*vhe);
}

/* Circular orbit velocity at a given altitude */
double vcirc(struct planet p, double r)
{
	return sqrt(p.gp / (p.r + r));
}

bool near(double a, double b)
{
	if (b == 0)
		return a == 0;
	return fabs(1.0 - a / b) < 1e-8;
}

static inline double fy(double x, double beta)
{
	return sqrt(1 - beta * beta * (1 - x * x)) * (beta < 0 ? -1.0 : 1.0);
}

static inline double arcoth(double x)
{
	return 0.5 * log((x + 1.0) / (x - 1.0));
}

int main(void)
{
	/* An example for testing: first Mars window after 1951.  Should be 6404m/s between 300km orbits */
	/* Depart vhe should be ~3952m/s, for ~3893m/s depart deltaV */
	double ut = 29033752-epoch;
	double tt = 298 * 86400;
	bool prograde = true;

	struct vector r1 = position_at_ut(planets[2], ut); // Earth depart
	struct vector r2 = position_at_ut(planets[3], ut + tt); // Mars arrive
	struct vector pv1 = velocity_at_ut(planets[2], ut); // Earth depart
	struct vector pv2 = velocity_at_ut(planets[3], ut + tt); // Mars arrive

	print_vector(r1);
	print_vector(pv1);
	print_vector(r2);
	print_vector(pv2);
	fprintf(stderr, "\n");

	/* Lambert solver, translated from https://github.com/alexmoon/ksp/blob/gh-pages/src/lambert.coffee */
	/* For simplicity, we use linear bisection where the original used Brent's Method */
	struct vector d = difference(r2, r1);
	double c = magnitude(d);
	double mr1 = magnitude(r1), mr2 = magnitude(r2);
	double m = mr1 + mr2 + c, n = mr1 + mr2 - c;
	double alpha = acos(dot(r1, r2) / (mr1 * mr2));
	double w = wedge(r1, r2).z;
	if ((w < 0) == prograde)
		alpha = M_PI * 2.0 - alpha;
	double beta = sqrt(n / m);
	if (alpha > M_PI)
		beta = -beta;
	double tn = 4.0 * tt * sqrt(sungp / pow(m, 3));
	double tnp = 2.0 / (3.0 * (1.0 - pow(beta, 3)));

	double x, y;

	if (near(tn, tnp)) { /* Parabolic */
		x = 1.0;
		y = beta < 0.0 ? -1.0 : 1.0;
	} else if (tn < tnp) { /* Hyperbolic */
		/* search for upper bound */
		double left, right = 1.0, val;
		do {
			left = right;
			x = right = left * 2.0;
			y = fy(x, beta);
			double g = sqrt(x * x - 1), h = sqrt(y * y - 1);
			val = (-arcoth(x/g) +arcoth(h/y) + x*g - y*h) / pow(g, 3) - tn;
		} while (val >= 0);
		if isnan(val) { /* There is apparently a bug that causes this to happen */
			fprintf(stderr, "NaN detected!  %g %g\n", left, right);
			return 1;
		}
		while (fabs(right - left) > 1e-8) {
			x = (left + right) / 2.0;
			y = fy(x, beta);
			double g = sqrt(x * x - 1), h = sqrt(y * y - 1);
			double val = (-arcoth(x/g) +arcoth(h/y) + x*g - y*h) / pow(g, 3) - tn;
			if (val < 0)
				left = x;
			else
				right = x;
		}
	} else { /* Elliptical */
		double tnme = acos(beta) + beta * sqrt(1 - pow(beta, 2));
		if (near(tn, tnme)) {
			x = 0.0;
			y = sqrt(1.0 - pow(beta, 2));
			if (beta < 0)
				y = -y;
		} else {
			double left, right;
			if (tn > tnme) {
				left = -1.0 + 1e-10;
				right = 0;
			} else {
				left = 0;
				right = 1.0 - 1e-10;
			}
			double yl = fy(left, beta),
			       yr = fy(right, beta),
			       gl = sqrt(1 - left * left),
			       gr = sqrt(1 - right * right),
			       hl = sqrt(1 - yl * yl),
			       hr = sqrt(1 - yr * yr),
			       vl = (M_PI / 2.0 - atan(left / gl) - atan(hl / yl) - left * gl + yl * hl) / pow(gl, 3) - tn,
			       vr = (M_PI / 2.0 - atan(right / gr) - atan(hr / yr) - right * gr + yr * hr) / pow(gr, 3) - tn;
			if ((vl < 0) == (vr < 0)) { /* can't happen */
				fprintf(stderr, "Solution out of range %g %g\n", vl, vr);
				return 1;
			} else if (vr < 0) { /* Swap endpoints */
				double tmp = left;
				left = right;
				right = tmp;
			}
			while (fabs(right - left) > 1e-10) {
				x = (left + right) / 2.0;
				y = fy(x, beta);
				double g = sqrt(1 - x * x), h = sqrt(1 - y * y);
				double val = (M_PI / 2.0 - atan(x / g) - atan(h / y) - x * g + y * h) / pow(g, 3) - tn;
				if (val < 0)
					left = x;
				else
					right = x;
			}
		}
	}


	double vc = sqrt(sungp) * (y / sqrt(n) + x / sqrt(m));
	double vr = sqrt(sungp) * (y / sqrt(n) - x / sqrt(m));
	struct vector ec = scale(d, vc / c);
	struct vector v1 = sum(ec, scale(r1, vr / mr1));
	struct vector v2 = difference(ec, scale(r2, vr / mr2));

	/* Heliocentric velocities */
	print_vector(v1);
	print_vector(v2);

	/* Relative velocities of escape and arrival */
	struct vector rv1 = difference(v1, pv1);
	struct vector rv2 = difference(v2, pv2);
	printf("Depart: %g\n", magnitude(rv1));
	print_vector(rv1);
	printf("Arrive: %g\n", magnitude(rv2));
	print_vector(rv2);

	/* Relative velocities at periapsis */
	double vp1 = v_from_vhe(planets[2], magnitude(rv1), 300e3);
	double vp2 = v_from_vhe(planets[3], magnitude(rv2), 300e3);
	printf("periV:\n");
	printf("  dep: %g\n", vp1);
	printf("  arr: %g\n", vp2);

	/* deltas V assuming optimal inclinations */
	double dv1 = vp1 - vcirc(planets[2], 300e3);
	double dv2 = vp2 - vcirc(planets[3], 300e3);
	printf("deltaV:\n");
	printf("   dep: %g\n", dv1);
	printf("   arr: %g\n", dv2);
	printf(" total: %g\n", dv1 + dv2);
	return 0;
}
