#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include "../utils/Random.hh"
#include "../utils/utils.hh"
#include "fft_mkl.hh"
using namespace std;


/*
 * 1D complex array
 * simple fn with two harmonics
 * look at coeffs
 */
void testc2c_1(int n){
	double *space;
	space = (double *)MKL_malloc(sizeof(double)*n*2, 32);
	double *v=space;
	for(int i=0; i < n; i++){
		double x = 2*PI*i*1.0/n;
		v[2*i] = 7*cos(4*x);
		v[2*i+1] = -3*cos(x);
	}
	
	FFTc2c fft(n, 1, n, 1);
	fft.fwd(v);
	
	cout<<fixed;
	cout<<setw(25)<<"MKL c2c"<<" easy fn = 7*cos(4x)-3i*cos(x)"<<endl;
	cout<<endl<<setw(10)<<"i"<<setw(20)<<"Fourier coeff"<<endl;
	for(int i=0; i < n; i++){
		cout<<setw(10)<<i<<setw(5)<<"("<<v[2*i]<<"+i*"<<v[2*i+1]
		    <<")"<<endl;
	}
	MKL_free(space);
}

/*
 * 1D complex array
 * random entries
 * compare fwd+bwd with original array
 */
void testc2c_2(int n){
	double *space;
	space = (double *)MKL_malloc(sizeof(double)*n*4, 32);
	double *v = space;
	double *w = v+2*n;
	Random rng;
	for(int i=0; i < n; i++){
		v[2*i] = rng.randn();
		v[2*i+1] = rng.randn();
		w[2*i] = v[2*i];
		w[2*i+1] = v[2*i+1];
		
	}

	cout<<setw(50)<<"c2c FFT test with random numbers"<<endl;
	cout<<setw(50)<<"n: "<<n<<endl;
	cout<<setw(25)<<"inf norm of v: "<<array_max(v, 2*n)<<endl;
	FFTc2c fft(n, 1, n, 1);
	fft.fwd(v);
	fft.bwd(v);
	array_diff(v, w, 2*n);
	cout<<setw(25)<<"n: "<<setw(10)<<n<<endl;
	cout<<setw(25)<<"error: "<<setw(10)
	    <<array_max(v, 2*n)/array_max(w, 2*n)
	    <<endl;
		
	MKL_free(space);
}

/*
 * several 1D arrays
 * blocked in sequence
 * fwd + bwd compared to original array
 */
void testc2c_3(int n, int count){
	double *space;
	space = (double *)MKL_malloc(sizeof(double)*n*count*4, 32);
	double *v = space;
	double *w = v+2*n*count;
	Random rng;
	for(int i=0; i < n*count; i++){
		v[2*i] = rng.randn();
		v[2*i+1] = rng.randn();
		w[2*i] = v[2*i];
		w[2*i+1] = v[2*i+1];
	}

	cout<<setw(50)<<"c2c FFT test with random entries"<<endl;
	cout<<setw(50)<<"n: "<<n<<endl;
	cout<<setw(50)<<"count of blocks: "<<count<<endl;
	cout<<setw(25)<<"inf norm of v: "<<array_max(v, 2*n*count)<<endl;
	FFTc2c fft(n, count, n, 1);
	fft.fwd(v);
	fft.bwd(v);
	array_diff(v, w, 2*n*count);
	cout<<setw(25)<<"n: "<<setw(10)<<n<<endl;
	cout<<setw(25)<<"error: "<<setw(10)
	    <<array_max(v, 2*n*count)/array_max(w, 2*n*count)
	    <<endl;
		
	MKL_free(space);
}

/*
 * interleaved 1D arrays
 * lda >= count must hold
 * superposed harmonics in each vector
 */
void testc2c_4(int n, int count, int lda){
	assrt(lda>=count);
	double *space;
	space = (double *)MKL_malloc(sizeof(double)*n*lda*4, 32);
	double *v = space;
	double *w = v+2*n*lda;
	for(int i=0; i < n*lda; i++){
		if((i%lda)>2*count)
			v[2*i] = v[2*i+1] = -1;
		else{
			double j = i/lda;
			double x = 2.0*PI*j/n;
			int blk = i%lda;
			v[2*i] = 2.0*(blk+1)*cos(x);
			v[2*i+1] = -4.0*(blk+1)*sin(3*x);
		}
		w[2*i] = v[2*i];
		w[2*i+1] = v[2*i+1];
	}
	
	cout<<setw(50)<<"c2c FFT test with superposed harmonics"<<endl;
	cout<<setw(50)<<"n: "<<n<<endl;
	cout<<setw(50)<<"number of interleaved blocks: "<<count<<endl;
	cout<<setw(50)<<"lda: "<<lda<<endl;
	cout<<setw(50)<<"block b: "<<"2(b+1)*cos(x)-4i*(b+1)*sin(3x)"<<endl;
	cout<<setw(25)<<"inf norm of v: "<<array_max(v, 2*n*count)<<endl;
	FFTc2c fft(n, count, 1, lda);
	fft.fwd(v);
	cout<<fixed;
	int blk = count/2;
	cout<<setw(25)<<"blk: "<<blk<<endl;
	for(int i=0; i < 8*lda; i+=lda)
		cout<<setw(25)<<"Fourier coeff blk : "
		    <<"("<<v[2*i+2*blk]<<","<<v[2*i+2*blk+1]<<")"<<endl;
	fft.bwd(v);
	array_diff(v, w, 2*n*count);
	cout<<setw(25)<<"n: "<<setw(10)<<n<<endl;
	cout<<setw(25)<<"error: "<<setw(10)<<scientific
	    <<array_max(v, 2*n*count)/array_max(w, 2*n*count)
	    <<endl;

	MKL_free(space);
}

int main(){
	testc2c_4(1024, 256, 261);
}
