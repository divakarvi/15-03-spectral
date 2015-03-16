#ifndef __BVPSolvei13July2013__
#define __BVPSolvei13July2013__
#include "../banded/tridiag.hh"
#include "../fft/trig.hh"

class BVPSolvei{
private:
	int state; //initially 0
	double a;
	int M, Me, Mo;//Mo is "M odd" not "M outer" in this class
	TriSolve even_tri;
	TriSolve odd_tri;
	const double *cheb;
	double *space, *h1, *h2, *d1, *d2;
	DCT *trig;
public:
	BVPSolvei(double aa, int MM); //state moves to 1
	~BVPSolvei();
	/*
	 * initcheb(), initdct(), inithmg() must be called just once
	 * and in that order
	 */
	void initcheb(double *chebi); //state moves to 2
	void initdct(DCT& dcti); //state moves to 3
	void inithmg(); //state move up to 4
	/*
	 * particular soln(s) of (D^2 - a^2)u = f
	 */
	void solvep(const double *restrict f,
		    double *restrict u, double *restrict du, int nrhs);
	/*
	 * particular soln(s) of (D^2 - a^2)u = f + dg/dy
	 */
	void solvepx(const double *restrict f, double *restrict g,
		     double *restrict u, double *restrict du, int nrhs);
	double *geth1(){ return h1;}
	double *geth2(){ return h2;}
	double *getd1(){ return d1;}
	double *getd2(){ return d2;}
};

#endif

