#ifndef __BVP22Mar15DV__
#define __BVP22Mar15DV__
#include "../banded/banded.hh"
#include "../fft/trig.hh"
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

class BVP4si{
private:
	double a, b; //parameters in (D^4 + a D^2 + b)u = f
	int M, Me, Mo;
	BandedSolve even;
	BandedSolve odd;
	DCT dct;
	double *space;
public:
	BVP4si(double aa, double bb, int MM);
	~BVP4si();
	/*
	 * f[0...M]
	 * u[0...M] is particular soln of (D^4 + a D^2 + b)u = f
	 * with T0, T1, T2, T3 cheb coeffs equal to zero
	 * dup1 and dum1 are the values of Du at +1 and -1, respectively
	 */
	void solvep(const double *restrict f, double *restrict u,
		    double& dup1, double& dum1);
};

#endif
