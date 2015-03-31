#ifndef __BVP22Mar15DV__
#define __BVP22Mar15DV__
#include "../banded/banded.hh"
#include "../fft/trig.hh"
#include "../legendre/legendre.hh"
#include "bvpi.hh"

/*
 * (D^2-a^2)*(D^2-b^2)u = f
 */
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

enum DCTONOFF {DCTON, DCTOFF};

/*
 * (D^4 + a D^2 + b)u = f
 */
class BVP4si{
private:
	double a, b; //parameters in (D^4 + a D^2 + b)u = f
	int M, Me, Mo;
	BandedSolve even;
	BandedSolve odd;
	DCT dct;
	double *space;
	double *h1, *h2, *h3, *h4;
	double H[16];
	int ipiv[4];
public:
	BVP4si(double aa, double bb, int MM);
	~BVP4si();
	/*
	 * f[0...M]
	 * u[0...M] is particular soln of (D^4 + a D^2 + b)u = f
	 * with T0, T1, T2, T3 cheb coeffs equal to zero
	 * dup1 and dum1 are the values of Du at +1 and -1, respectively
	 * if flag = DCTOFF a final bwd DCT is omitted and the u returned
	 * is the Cheby series
	 */
	void solvep(const double *restrict f, double *restrict u,
		    double& dup1, double& dum1, enum DCTONOFF flag = DCTON);
	void solve(const double *restrict f, double *restrict u);
};

/*
 * (D^4 - b*D^2 + a)u = f 
 */
class BVP4pg{
private:
	double a,b;
	int M, Me, Mo;
	CholeskyBanded even;
	CholeskyBanded odd;
	double *space;
	FLTrans flt;
public:
	BVP4pg(double aa, double bb, int MM);
	~BVP4pg();
	void solve(const double *restrict f, double *restrict u);
};

#endif
