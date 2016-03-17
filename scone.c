#include <stdio.h>
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

int main(void)
{
	/* An example for testing */
	double ut = 1661593873.0 - epoch; // closest approach in modern times
	struct vector e = position_at_ut(planets[2], ut);
	struct vector m = position_at_ut(planets[3], ut);

	printf("%g, %g, %g\n", e.x, e.y, e.z);
	printf("%g, %g, %g\n", m.x, m.y, m.z);
	printf("%g\n", magnitude(difference(e, m)));
	return 0;
}
