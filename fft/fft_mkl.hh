#ifndef __MKLFourier24June2013__
#define __MKLFourier24June2013__
#include <mkl.h>
#include <mkl_dfti.h>
class FFTc2c{
private:
	DFTI_DESCRIPTOR_HANDLE handle;
public:
	/*
	 * n = size of transform (unit=complex)
	 * count = number of vectors to be transformed
	 * dist = dist betweent first entries of succ vectors (unit=complex)
	 * stride = dist between two entries of same vector (unit=complex)
	 */
	FFTc2c(long int n, int count, int dist, int stride); 
	~FFTc2c();
	//size>=2*(dist*(count-1)+stride*n)
	void fwd(double *v){
		DftiComputeForward(handle, v);
	};
	//size>=2*(dist*(count-1)+stride*n)
	void bwd(double *vf){
		DftiComputeBackward(handle, vf);
	};
};

class FFTr2c{
private:
	DFTI_DESCRIPTOR_HANDLE fhandle, bhandle;
public:
	/*
	 * n = size of transform (unit=double)
	 * count = number of transforms
	 * the vectors or blocks must be sequentially arranged
	 * interleaving is not allowed
	 * dist = distance between first entries of succ vectors (unit=double)
	 * dist >= 2((n/2)+1) is required
	 */
	FFTr2c(int n, int count, int dist); 
	~FFTr2c();
	//size >= count * dist
	void fwd(double *v){
		DftiComputeForward(fhandle, v);
	};
	//size >= count * dist
	void bwd(double *vf){
		DftiComputeBackward(bhandle, vf);
	};
};

class FFT2D{
private:
	DFTI_DESCRIPTOR_HANDLE fhandle;
	DFTI_DESCRIPTOR_HANDLE bhandle;
public:
	/*
	 * m = num of rows  
	 * n = num of cols (real data) (row major)
	 * lda = m x n is stuffed inside m x lda array
	 * lda >= 2*(n/2+1)
	 * lda must be even
	 * count = num of transforms
	 */
	FFT2D(int m, int n, int lda, int count);
	~FFT2D();
	//size >= m*lda*count 
	void fwd(double *v){
		DftiComputeForward(fhandle, v);
	}
	//size >= m*lda*count 
	void bwd(double *vf){
		DftiComputeBackward(bhandle, vf);
	}
};

#endif
