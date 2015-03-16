#include "fft2D_thrd.hh"
#include "../utils/utils.hh"
#include <iostream>
using namespace std;

/*
 * 0=fst[0]<=...<=fst[P]=n at output
 * the block owned by p is fst[p]<=i < fst[p+1]
 * number of items owned by p is fst[p+1]-fst[p]
 * blocks are as nearly equal as possible
 * P = nprocs
 */
static void BlockDivide(long n, int P, long *fst){
	long  Q = n/P;
	long  R = n-Q*P;
	fst[0] = 0;
	for(int p=0; p < R; p++)
		fst[p+1] = fst[p] + (Q+1);
	for(int p=R; p < P; p++)
		fst[p+1] = fst[p]+Q;
}

FFT2D_thrd::FFT2D_thrd(int Ti, int mi, int ni, int ldai, int counti)
:T(Ti), m(mi), n(ni), lda(ldai), count(counti)
{
	assrt(m>0);
	assrt(n>0);
	assrt(lda>=2*(n/2+1));
	assrt(lda%2==0);
	assrt(count>0);
	fstmc = new long[T+1];
	fstnby2p1 = new long[T+1];
	typedef FFTr2c *fftr2ctemp;
	r2clist = new fftr2ctemp[T+1];
	typedef FFTc2c *fftc2ctemp;
	c2clist = new fftc2ctemp[T+1];
	
	BlockDivide(m*count, T, fstmc);
	BlockDivide(n/2+1, T, fstnby2p1);

#pragma omp parallel				\
	num_threads(T)				\
	default(shared)
	{
		int t = omp_get_thread_num();
		r2clist[t] = new FFTr2c(n, fstmc[t+1]-fstmc[t], lda);
		c2clist[t] = new FFTc2c(m, fstnby2p1[t+1]-fstnby2p1[t], 
					1, lda/2);
	}
}

FFT2D_thrd::~FFT2D_thrd(){
	for(int t=0; t < T; t++){
		delete r2clist[t];
		delete c2clist[t];
	}
	delete[] c2clist;
	delete[] r2clist;
	delete[] fstnby2p1;
	delete[] fstmc;
}

void FFT2D_thrd::fwd(double *v){
#pragma omp parallel				\
	num_threads(T)				\
	default(shared)
	{
		int t = omp_get_thread_num();
		double *vt = v + fstmc[t]*lda;
		r2clist[t]->fwd(vt);
	}
#pragma omp parallel				\
	num_threads(T)				\
	default(shared)				
	{
		int t = omp_get_thread_num();
		double *vt = v + fstnby2p1[t]*2;
		for(int c=0; c < count; c++){
			c2clist[t]->fwd(vt);
			vt += m*lda;
		}
			
	}
}

void FFT2D_thrd::bwd(double *v){
#pragma omp parallel				\
	num_threads(T)				\
	default(shared)
	{
		int t = omp_get_thread_num();
		double *vt = v + fstnby2p1[t]*2;
		for(int c=0; c < count; c++){
			c2clist[t]->bwd(vt);
			vt += m*lda;
		}
			
	}
#pragma omp parallel				\
	num_threads(T)				\
	default(shared)
	{
		int t = omp_get_thread_num();
		double *vt = v + fstmc[t]*lda;
		r2clist[t]->bwd(vt);
	}
}
