#include <math.h>

struct vector {
	double x;
	double y;
	double z;
};

struct vector difference(struct vector v, struct vector w);
double magnitude(struct vector v);

/* All angles in radians */
struct vector rotateXY(struct vector v, double theta);
struct vector rotateXZ(struct vector v, double theta);
struct vector rotateYZ(struct vector v, double theta);
