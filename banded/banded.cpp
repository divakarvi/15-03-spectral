#include "banded.hh"
#include "../utils/utils.hh"
#include <mkl.h>

BandedSolve::BandedSolve(int ni)
	:n(ni), state(0){
	space = (double *)MKL_malloc(sizeof(double)*7*n
					     + sizeof(int)*n, 16);
	A = space;
	ipiv = (int *)(A+7*n);
}

BandedSolve::~BandedSolve(){
	MKL_free(space);
}

void BandedSolve::setu2(double *u2){
	assrt(state==0);
	state = 1;
	int i = 2;
	int lda = 7;
	for(int j = 2; j < n; j++)
		A[i+j*lda] = u2[j-2];
}


void BandedSolve::setu1(double *u1){
	assrt(state==1);
	state = 2;
	int i = 3;
	int lda = 7;
	for(int j = 1; j < n; j++)
		A[i+j*lda] = u1[j-1];
}

void BandedSolve::setdiag(double *d){
	assrt(state==2);
	state = 3;
	int i = 4;
	int lda = 7;
	for(int j = 0; j < n; j++)
		A[i+j*lda] = d[j];
}

void BandedSolve::setl1(double *l1){
	assrt(state==3);
	state = 4;
	int i = 5;
	int lda = 7;
	for(int j = 0; j < n-1; j++)
		A[i+j*lda] = l1[j];
}

void BandedSolve::setl2(double *l2){
	assrt(state==4);
	state = 5;
	int i = 6;
	int lda = 7;
	for(int j = 0; j < n-2; j++)
		A[i+j*lda] = l2[j];
}

void BandedSolve::factor(){
	assrt(state==5);
	state = 6;
	int klu = 2;
	int lda = 7;
	int info;
	dgbtrf_(&n, &n, &klu, &klu, A, &lda, ipiv, &info);
	assrt(info==0);
}

void BandedSolve::solve(double *f, int nrhs){
	assrt(state==6);
	assrt(nrhs > 0);
	char trans[3] = "N";
	int klu = 2;
	int lda = 7;
	int ldf = n;
	int info;
	dgbtrs_(trans, &n, &klu, &klu, &nrhs, A, &lda, ipiv, f, &ldf, &info);
	assrt(info==0);
}
