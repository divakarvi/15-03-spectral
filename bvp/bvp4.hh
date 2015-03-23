#ifndef __BVP22Mar15DV__
#define __BVP22Mar15DV__
#include "../trig/trig.hh"
#include "bvpi.hh"

class BVP4fac{
private:
	int M;
	BVPSolvei alpha;
	BVPSolvei beta;

	double *chebi;
	DCT dct;

	double *h1, *h2, *h3, *h4, *du, *w;
	double H[16];
	int ipiv[4];
public:
	BVP4fac(double a, double b,  int MM);
	~BVP4fac();
	void solve(const double *restrict f, double *restrict u);

};

#endif
