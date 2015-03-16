#include "fft_mkl.hh"
#include "../utils/utils.hh"

FFTc2c::FFTc2c(long int n, int count, int dist, int stride)
{
	assrt(n > 0);
	assrt(count > 0);
	
	DftiCreateDescriptor(&handle, DFTI_DOUBLE, DFTI_COMPLEX, 1, n);
	DftiSetValue(handle, DFTI_FORWARD_SCALE, 1.0/n);
	DftiSetValue(handle, DFTI_PLACEMENT, DFTI_INPLACE);
	DftiSetValue(handle, DFTI_NUMBER_OF_TRANSFORMS,count);
	DftiSetValue(handle, DFTI_INPUT_DISTANCE, dist);
	DftiSetValue(handle, DFTI_OUTPUT_DISTANCE, dist);
	MKL_LONG str[2];
	str[0] = 0;
	str[1] = stride;
	DftiSetValue(handle, DFTI_INPUT_STRIDES, str);
	DftiSetValue(handle, DFTI_OUTPUT_STRIDES, str);
	DftiCommitDescriptor(handle);
}

FFTc2c::~FFTc2c(){
	DftiFreeDescriptor(&handle);
}

FFTr2c::FFTr2c(int n, int count, int dist){
	assrt(n > 0);
	assrt(count > 0);
	assrt(dist%2==0);
	assrt(dist >= 2*(n/2+1));
		
	/* 
	 * strides are set to 1 by default
	 */
	DftiCreateDescriptor(&fhandle,DFTI_DOUBLE,DFTI_REAL,1,n);
	DftiSetValue(fhandle, DFTI_PACKED_FORMAT, DFTI_CCE_FORMAT);
	DftiSetValue(fhandle, DFTI_CONJUGATE_EVEN_STORAGE, 
		      DFTI_COMPLEX_COMPLEX);
	DftiSetValue(fhandle, DFTI_FORWARD_SCALE, 1.0/n);
	DftiSetValue(fhandle, DFTI_PLACEMENT, DFTI_INPLACE);
	DftiSetValue(fhandle, DFTI_NUMBER_OF_TRANSFORMS,count);
	DftiSetValue(fhandle, DFTI_INPUT_DISTANCE, dist);
	DftiSetValue(fhandle, DFTI_OUTPUT_DISTANCE, dist/2);
	DftiCommitDescriptor(fhandle);

	DftiCreateDescriptor(&bhandle,DFTI_DOUBLE,DFTI_REAL,1,n);
	DftiSetValue(bhandle, DFTI_PACKED_FORMAT, DFTI_CCE_FORMAT);
	DftiSetValue(bhandle, DFTI_CONJUGATE_EVEN_STORAGE,
		     DFTI_COMPLEX_COMPLEX);
	DftiSetValue(bhandle, DFTI_PLACEMENT, DFTI_INPLACE);
	DftiSetValue(bhandle, DFTI_NUMBER_OF_TRANSFORMS,count);
	DftiSetValue(bhandle, DFTI_INPUT_DISTANCE, dist/2);
	DftiSetValue(bhandle, DFTI_OUTPUT_DISTANCE, dist);
	DftiCommitDescriptor(bhandle);
}

FFTr2c::~FFTr2c(){
	DftiFreeDescriptor(&fhandle);
	DftiFreeDescriptor(&bhandle);
}

FFT2D::FFT2D(int m, int n, int lda, int count){
	assrt(m > 0);
	assrt(n > 0);
	assrt(lda >= 2*(n/2+1));
	assrt(lda%2 == 0);
	assrt(count > 0);
	
	long int nlist[2];
	nlist[0] = m;
	nlist[1] = n;
	DftiCreateDescriptor(&fhandle, DFTI_DOUBLE, DFTI_REAL, 2, nlist);
	DftiSetValue(fhandle, DFTI_PACKED_FORMAT, DFTI_CCE_FORMAT);
	DftiSetValue(fhandle, DFTI_CONJUGATE_EVEN_STORAGE, 
		     DFTI_COMPLEX_COMPLEX);
	DftiSetValue(fhandle, DFTI_FORWARD_SCALE, 1.0/m/n);
	DftiSetValue(fhandle, DFTI_PLACEMENT, DFTI_INPLACE);
	DftiSetValue(fhandle, DFTI_NUMBER_OF_TRANSFORMS,count);
	int dist = m*lda;
	DftiSetValue(fhandle, DFTI_INPUT_DISTANCE, dist);
	DftiSetValue(fhandle, DFTI_OUTPUT_DISTANCE, dist/2);
	MKL_LONG stride[3];
	stride[0] = 0;
	stride[1] = lda;
	stride[2] = 1;
	DftiSetValue(fhandle, DFTI_INPUT_STRIDES, stride);
	stride[0] = 0;
	stride[1] = lda/2;
	stride[2] = 1;
	DftiSetValue(fhandle, DFTI_OUTPUT_STRIDES, stride);
	DftiCommitDescriptor(fhandle);

	DftiCreateDescriptor(&bhandle, DFTI_DOUBLE, DFTI_REAL, 2, nlist);
	DftiSetValue(bhandle, DFTI_PACKED_FORMAT, DFTI_CCE_FORMAT);
	DftiSetValue(bhandle, DFTI_CONJUGATE_EVEN_STORAGE, 
		     DFTI_COMPLEX_COMPLEX);
	DftiSetValue(bhandle, DFTI_PLACEMENT, DFTI_INPLACE);
	DftiSetValue(bhandle, DFTI_NUMBER_OF_TRANSFORMS,count);
	dist = m*lda;
	DftiSetValue(bhandle, DFTI_INPUT_DISTANCE, dist/2);
	DftiSetValue(bhandle, DFTI_OUTPUT_DISTANCE, dist);
	stride[0] = 0;
	stride[1] = lda/2;
	stride[2] = 1;
	DftiSetValue(bhandle, DFTI_INPUT_STRIDES, stride);
	stride[0] = 0;
	stride[1] = lda;
	stride[2] = 1;
	DftiSetValue(bhandle, DFTI_OUTPUT_STRIDES, stride);
	DftiCommitDescriptor(bhandle);
}

FFT2D::~FFT2D(){
	DftiFreeDescriptor(&fhandle);
	DftiFreeDescriptor(&bhandle);
}
