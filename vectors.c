#include "vectors.h"

#include <math.h>

struct vector difference(struct vector v, struct vector w)
{
	return (struct vector) {
		.x = v.x - w.x,
		.y = v.y - w.y,
		.z = v.z - w.z,
	};
}

double magnitude(struct vector v)
{
	return sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
}

struct vector rotateXY(struct vector v, double theta)
{
	struct vector w;
	w.x = v.x * cos(theta) - v.y * sin(theta);
	w.y = v.x * sin(theta) + v.y * cos(theta);
	w.z = v.z;
	return w;
}

struct vector rotateXZ(struct vector v, double theta)
{
	struct vector w;
	w.x = v.x * cos(theta) - v.z * sin(theta);
	w.y = v.y;
	w.z = v.x * sin(theta) + v.z * cos(theta);
	return w;
}

struct vector rotateYZ(struct vector v, double theta)
{
	struct vector w;
	w.x = v.x;
	w.y = v.y * cos(theta) - v.z * sin(theta);
	w.z = v.y * sin(theta) + v.z * cos(theta);
	return w;
}
