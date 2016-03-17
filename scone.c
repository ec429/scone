#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include "planets.h"
#include "vectors.h"

double ea_from_ma(double ma, double ecc)
{
	double left = 0, right = M_PI * 2.0, mid, val;

	ma = fmod(ma, 2 * M_PI);
	while (right - left > 1e-8) {
		mid = (left + right) / 2.0;
		val = mid - ecc * sin(mid);
		if (val < ma)
			left = mid;
		else
			right = mid;
	}
	return mid;
}

struct vector position_at_ut(struct planet p, double ut)
{
	double period = 2.0 * M_PI * sqrt(pow(p.sma, 3) / sungp);
	double mean_anomaly = 2.0 * M_PI * (ut - epoch) / period;
	double ecc_anomaly = ea_from_ma(mean_anomaly, p.ecc);
	double true_anomaly = 2.0 * atan2(sqrt(1.0 + p.ecc) * sin(ecc_anomaly / 2.0),
					  sqrt(1.0 - p.ecc) * cos(ecc_anomaly / 2.0));
	struct vector v = {p.sma / (1.0 + p.ecc * cos(true_anomaly)), 0, 0};

	v = rotateXY(v, true_anomaly + p.ape);
	v = rotateYZ(v, p.inc);
	v = rotateXY(v, p.lan);
	return v;
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
	/* An example for testing: 300-day transfer starting at game start */
	double ut = -epoch; // 19510101T0000Z
	double tt = 300 * 86400; // 300 days
	bool prograde = true;
	struct vector r1 = position_at_ut(planets[2], ut); // Earth depart
	struct vector r2 = position_at_ut(planets[3], ut + tt); // Mars arrive

	print_vector(r1);
	print_vector(r2);

	/* Lambert solver, translated from https://github.com/alexmoon/ksp/blob/gh-pages/src/lambert.coffee */
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
			val = (-arcoth(x/g) +arcoth(h/y) - x*g + y*h) / pow(g, 3) - tn;
		} while (val >= 0);
		if isnan(val) { /* There is apparently a bug that causes this to happen */
			fprintf(stderr, "NaN detected!  %g %g\n", left, right);
			return 1;
		}
		while (fabs(right - left) > 1e-8) {
			x = (left + right) / 2.0;
			y = fy(x, beta);
			double g = sqrt(x * x - 1), h = sqrt(y * y - 1);
			double val = (-arcoth(x/g) +arcoth(h/y) - x*g + y*h) / pow(g, 3) - tn;
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
				left = -1.0 + 1e-8;
				right = 0;
			} else {
				left = 0;
				right = 1.0 - 1e-8;
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
			while (fabs(right - left) > 1e-8) {
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


	double vc = sqrt(sungp) * (y / sqrt(n) + 1.0 / sqrt(m));
	double vr = sqrt(sungp) * (y / sqrt(n) - 1.0 / sqrt(m));
	struct vector ec = scale(d, vc / c);
	struct vector v1 = sum(ec, scale(r1, vr / mr1));
	struct vector v2 = sum(ec, scale(r2, vr / mr2));

	/* Heliocentric velocities */
	print_vector(v1);
	print_vector(v2);

	return 0;
}
