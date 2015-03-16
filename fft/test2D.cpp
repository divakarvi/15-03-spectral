#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include "../utils/Random.hh"
#include "../utils/utils.hh"
#include "fft_mkl.hh"
using namespace std;

/*
 * single set of 2D data
 * superposed harmonics
 * check Fourier coeffs
 */
void test2D_1(int m, int n){
	assrt(m > 10);
	assrt(n > 10);
	int lda = 2*(n/2+1);
	double *space;
	space = (double *)MKL_malloc(sizeof(double)*2*m*lda, 32);
	double *v = space;
	double *w = v + m*lda;
	for(int i=0; i < m; i++){
		for(int j=0; j < n; j++){
			double y = 2.0*PI*i/m;
			double x = 2.0*PI*j/n;
			v[i*lda+j] = 2.0*cos(x+y)-4.0*sin(3*x-2*y);
			w[i*lda+j] = v[i*lda+j];
		}
		for(int j=n; j < lda; j++){
			v[i*lda+j] = 1e100;
		}
	}
	
	cout<<setw(50)<<"2D FFT with count 1"<<endl;
	cout<<setw(50)<<"m: "<<m<<endl;
	cout<<setw(50)<<"n: "<<n<<endl;
	cout<<setw(50)<<"fn: "<<"2cos(x+y)-4sin(3x-2y)"<<endl;
	FFT2D fft(m, n, lda, 1);
	fft.fwd(v);
	for(int i=0; i < 4; i++)
		for(int j=0; j < 4; j++){
			cout<<setw(25)<<"ky: "<<i<<endl;
			cout<<setw(25)<<"kx: "<<j<<endl;
			cout<<fixed<<setw(25)<<"Fourier Coeff: "
			    <<"("<<v[i*lda+2*j]
			    <<","<<v[i*lda+2*j+1]<<")"<<endl;
			if(i==0)
				cout<<endl;
			else{
				cout<<setw(25)<<"ky: "<<-i<<endl;
				cout<<setw(25)<<"kx: "<<j<<endl;
				cout<<fixed<<setw(25)<<"Fourier Coeff: "
				    <<"("<<v[(m-i)*lda+2*j]
				    <<","<<v[(m-i)*lda+2*j+1]<<")"<<endl<<endl;
			}
		}
	cout<<endl<<endl;
	cout<<setw(50)<<"unused spots in CCE format not garbled"<<endl;
	//v[0*lda+1] = 1e100;
	//v[0*lda+n+1] = -1e100;
	cout<<setw(50)<<"garbling is dangerous with 2D FFTs"<<endl;
	
	fft.bwd(v);
	//fix unused spots before finding error
	for(int i=0; i < m; i++){
		v[i*lda+n] = w[i*lda+n] = 0;
		if(n%2==0)
			v[i*lda+n+1] = w[i*lda+n+1] = 0;
	}
	array_diff(v, w, m*lda);
	cout<<scientific<<setw(25)<<"error: "
	    <<array_max(v, m*lda)/array_max(w, m*lda)
	    <<endl<<endl;

	MKL_free(space);
}

/*
 * multiple sets of 2D data
 * superposed harmonics
 * check Fourier coeffs
 */
void test2D_2(int m, int n, int count){
	assrt(m > 10);
	assrt(n > 10);
	int lda = 2*(n/2+1)+10;//a little padding
	double *space;
	space = (double *)MKL_malloc(sizeof(double)*2*m*lda*count, 32);
	double *v = space;
	double *w = v + m*lda*count;
	for(int i=0; i < m*lda*count; i++)//for error check later
		v[i] = w[i] = 0;
	for(int c=0; c < count; c++)
		for(int i=0; i < m; i++){
			for(int j=0; j < n; j++){
				double y = 2.0*PI*i/m;
				double x = 2.0*PI*j/n;
				v[c*m*lda+i*lda+j] = 
					+(c+1)*2.0*cos(x+y)
					-(c+1)*4.0*sin(3*x-2*y);
				w[c*m*lda+i*lda+j] = v[c*m*lda+i*lda+j];
			}
			for(int j=n; j < lda; j++)
				v[j+i*lda+c*m*lda] = -1e100;
		}

	cout<<setw(50)<<"2D FFT with lda and count"<<endl;
	cout<<setw(50)<<"m: "<<m<<endl;
	cout<<setw(50)<<"n: "<<n<<endl;
	cout<<setw(50)<<"lda: "<<lda<<endl;
	cout<<setw(50)<<"blocks: "<<count<<endl;
	cout<<setw(50)<<"fn in blk b: "
	    <<"2(b+1)cos(x+y)-4(b+1)sin(3x-2y)"<<endl;
	FFT2D fft(m, n, lda, count);
	fft.fwd(v);
	int blk = count/2;
	cout<<setw(25)<<"block b: "<<blk<<endl;
	for(int i=0; i < 4; i++)
		for(int j=0; j < 4; j++){
			cout<<setw(25)<<"ky: "<<i<<endl;
			cout<<setw(25)<<"kx: "<<j<<endl;
			cout<<fixed<<setw(25)<<"Fourier Coeff: "
			    <<"("<<v[blk*m*lda+i*lda+2*j]
			    <<","<<v[blk*m*lda+i*lda+2*j+1]<<")"<<endl;
			if(i==0)
				cout<<endl;
			else{
				cout<<setw(25)<<"ky: "<<-i<<endl;
				cout<<setw(25)<<"kx: "<<j<<endl;
				cout<<fixed<<setw(25)<<"Fourier Coeff: "
				    <<"("<<v[blk*m*lda+(m-i)*lda+2*j]
				    <<","<<v[blk*m*lda+(m-i)*lda+2*j+1]
				    <<")"<<endl<<endl;
			}
		}
	cout<<endl<<endl;
	
	fft.bwd(v);

	//fix unused spots before finding error
	for(int c=0; c < count; c++)
		for(int i=0; i < m; i++)
			for(int j=n; j < lda; j++)
				v[c*m*lda+i*lda+j] = w[c*m*lda+i*lda+j] = 0;

	array_diff(v, w, count*m*lda);
	cout<<scientific<<setw(25)<<"error: "
	    <<array_max(v, count*m*lda)/array_max(w, count*m*lda)
	    <<endl<<endl;
	MKL_free(space);
}

int main(){
	test2D_1(1023, 2055);
	test2D_2(1023, 2059, 128);
}
