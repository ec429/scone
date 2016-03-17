#include <math.h>

struct vector {
	double x;
	double y;
	double z;
};

struct vector sum(struct vector v, struct vector w);
struct vector difference(struct vector v, struct vector w);
struct vector scale(struct vector v, double s);
double magnitude(struct vector v);

double dot(struct vector v, struct vector w);
struct vector wedge(struct vector v, struct vector w);

/* All angles in radians */
struct vector rotateXY(struct vector v, double theta);
struct vector rotateXZ(struct vector v, double theta);
struct vector rotateYZ(struct vector v, double theta);

void print_vector(struct vector v);
