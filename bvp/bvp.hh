#ifndef __BVPSolve16July2013__
#define __BVPSolve16July2013__
#include "bvpi.hh"
#include "../banded/banded.hh"
#include "../banded/tridiag.hh"
#include "../fft/trig.hh"

class BVPSolve{
private:
	int state; //initialized to zero
	/*
	 * a  = parameter in BVP (D*D-a*a)u = f
	 * Mi = M inner
	 * Mo = M outer
	 * tid = thread id (for claiming scr memory)
	 */
	double a;
	int Mi, Mo;
	BVPSolvei **bvpi;
	BandedSolve banded;
	double B[4]; //used if Mo=1
	const double *oy;
public:
	BVPSolve(double aa, int MMi, int MMo); //state moves to 1
	~BVPSolve();
	void init_oy(const double *oy); //state moves to 2
	void initcheb(double *chebi); //state moves to 3
	void initdct(DCT& dcti); //state moves to 4
	void initbanded(); //state moves to 5
	/*
	 * solve (D^2 - a^2)u = f + dg/dy
	 * f        = (Mi+1)*Mo doubles
	 * g        = (Mi+1)*Mo doubles (NULL implies g omitted)
	 * u        = soln
	 * du       = derv of soln
	 * plusone  = boundary value at +1
	 * minusone = boundary value at -1
	 * f may be overwritten by du
	 */
	void solve(const double *f, 
		   double *restrict u, double *du,
		   double plusone, double minusone,
		   double *restrict g = NULL);
	/*
	 * solve with 0 bndry conditions
	 */
	void solve(const double *f, double *restrict u, double *du){
		assrt(state==5);
		solve(f, u, du, 0.0, 0.0, NULL);
	}

	void solvex(const double *f, double *restrict g, 
		    double *restrict u, double *du){
		assrt(state==5);
		solve(f, u, du, 0.0, 0.0, g);
	}
	
	/*
	 * pa = homg soln satisfying u(1) = 1 u(-1) = 0;
	 * pb = homg soln satisfying u(1) = 0 u(-1) = 1;
	 * called from ks/ks.cpp
	 */
	void hmg10n01(double *restrict pa, double *restrict dpa, 
			 double *restrict pb, double *restrict dpb);
};

#endif
