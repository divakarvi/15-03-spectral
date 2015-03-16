#include "tridiag.hh"
#include "../utils/utils.hh"
#include <mkl.h>

TriSolve::TriSolve(int ni)
	:n(ni), state(0){
	space = (double *)MKL_malloc(sizeof(double)*
				     ((n-2)+(n-1)+n+(n-1)) //diagonals
				     + sizeof(int)*n, //ipiv
				     8);
	u2 = space;
	u1 =  u2 + (n-2);
	diag = u1 + (n-1);
	l1 = diag + n;
	ipiv = (int *)(l1+(n-1));
}

TriSolve::~TriSolve(){
	MKL_free(space);
}

void TriSolve::factor(){
	assrt(state==3);
	state = 4;
	int info;
	dgttrf_(&n, l1, diag, u1, u2, ipiv, &info);
	assrt(info==0);
}

void TriSolve::solve(double *f, int nrhs){
	assrt(state==4);
	assrt(nrhs>0);
	char trans[3] = "N";
	int ldf = n;
	int info;
	dgttrs_(trans, &n,  &nrhs, l1, diag, u1, u2, ipiv, f, &ldf, &info);
	assrt(info==0);
}
