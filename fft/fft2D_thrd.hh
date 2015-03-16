#ifndef __MKLFourierThrd28June2013__
#define __MKLFourierThrd28June2013__
#include "fft_mkl.hh"
#include <omp.h>

class FFT2D_thrd{
private:
	int T; //num of omp threads
	int m, n; //FFT dimensions
	int lda; //lda>=2*(n/2+1), lda even, unit = double
	int count; //num of FFTs
	long *fstmc; //array which divides m*count
	long *fstnby2p1; //array which divides (n/2+1)
	FFTr2c **r2clist;
	FFTc2c **c2clist;
public:
	/*
	 * Ti, mi, ni, ldai, counti ---> T, m, n, lda, count above
	 */
	FFT2D_thrd(int Ti, int mi, int ni, int ldai, int counti);
	~FFT2D_thrd();
	/*
	 * size of v >= m*lda*count (unit = double)
	 */
	void fwd(double *v);
	void bwd(double *v);
};
#endif
