#ifndef __TRIDIAG5JULY2013__
#define __TRIDIAG5JULY2013__
#include "../utils/utils.hh"

class TriSolve{//solve tridiagonal system
private:
	int n;
	int state;
	double *space;
	double *u1, *u2;
	double *diag;
	double *l1;
	int *ipiv;
public:
	TriSolve(int ni); //sets state to 0
	~TriSolve();
	/*
	 * the next four routines setu1(), setdiag(), setl1(), factor()
	 * must be called just once and in that order
	 */
	void setu1(double *restrict tmp){//size = (n-1) doubles 
		assrt(state==0);
		state = 1;
		for(int i=0; i < n-1; i++)
			u1[i] = tmp[i];
	}
	void setdiag(double *restrict tmp){//size = n doubles
		assrt(state==1);
		state = 2;
		for(int i=0; i < n; i++)
			diag[i] = tmp[i];
	}
	void setl1(double *restrict tmp){//size = (n-1) doubles
		assrt(state==2);
		state = 3;
		for(int i=0; i < n-1; i++)
			l1[i] = tmp[i];
	}
	void factor();
	/*
	 * calls to setu1(), setdiag(), setl1(), factor()
	 * must precede calls to solve()
	 * f    = (input) rhs at input
	 *      = (soln) soln at output
	 * nrhs = num of rhs (lda assumed to be n)
	 */
	void solve(double *f, int nrhs);
	double *getu1(){ return u1;}
	double *getdiag(){return diag;}
	double *getl1(){return l1;}
};
#endif
