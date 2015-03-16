#ifndef __BANDED5JULY2013__
#define __BANDED5JULY2013__
#include "../utils/utils.hh"

class BandedSolve{//solves banded matrices with 5 (2+1+2) diagonals
private:
	int n;
	double *space;//7xn matrix
	double *A;
	int *ipiv;
	int state;
public:
	BandedSolve(int ni);
	~BandedSolve();
	/*
	 * size of u2 = n-2 
	 * size of u1 = n-1 
	 * size of diag = n
	 * size of l1 = n-1
	 * size of l2 = n-2
	 * setu2(), setu1(), setdiag(), setl1(), setl2()
	 * must be called before factor() and in that order
	 */
	void setu2(double *u2);//size = (n-2) doubles
	void setu1(double *u1);//size = (n-1) doubles
	void setdiag(double *d);//size = n doubles
	void setl1(double *l1);//size = (n-1) doubles
	void setl2(double *l2);//size = (n-2) doubles
	void factor();
	/*
	 * f = (input) rhs
	 *   = (output) solution
	 * setu2(), setu1(), setdiag(), setl1(), setl2(), and factor()
	 * must be called in that order before calling solve()
	 */
	void solve(double *f, int nrhs);
	double *getA(){ return A;}
};
#endif
