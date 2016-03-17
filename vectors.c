#include "vectors.h"

#include <stdio.h>
#include <math.h>

struct vector sum(struct vector v, struct vector w)
{
	return (struct vector) {
		.x = v.x + w.x,
		.y = v.y + w.y,
		.z = v.z + w.z,
	};
}

struct vector difference(struct vector v, struct vector w)
{
	return (struct vector) {
		.x = v.x - w.x,
		.y = v.y - w.y,
		.z = v.z - w.z,
	};
}

struct vector scale(struct vector v, double s)
{
	return (struct vector) {
		.x = v.x * s,
		.y = v.y * s,
		.z = v.z * s,
	};
}

double magnitude(struct vector v)
{
	return sqrt(dot(v, v));
}

inline double dot(struct vector v, struct vector w)
{
	return v.x*w.x + v.y*w.y + v.z*w.z;
}

struct vector wedge(struct vector v, struct vector w)
{
	return (struct vector) {
		.x = v.y * w.z - v.z * w.y,
		.y = v.z * w.x - v.x * w.z,
		.z = v.x * w.y - v.y * w.x,
	};
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

void print_vector(struct vector v)
{
	printf("(%g, %g, %g)\n", v.x, v.y, v.z);
}
