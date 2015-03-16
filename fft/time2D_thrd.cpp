#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include "../utils/Random.hh"
#include "../utils/utils.hh"
#include "../utils/TimeStamp.hh"
#include "../utils/StatVector.hh"
#include "fft_mkl.hh"
#include "fft2D_thrd.hh"

/*
 * returns time for mxn fft using T threads
 */
double time2d(int m, int n,  int T){
	int pad = 0;
	int lda = 2*(n/2+1)+pad;
	int unitsize = m*lda*8;//in bytes
	int count = 8l*1000*1000*1000/unitsize;//8 GB
	double *space = (double *)MKL_malloc(1l*unitsize*count, 32);
	double *v = space;	
	FFT2D_thrd fft(T, m, n, lda, count);
	fft.fwd(v);//for numa
	for(long i=0; i < count*m*lda; i++)
		v[i] = (i-7.0)/(i*i+i+1);

	TimeStamp clk;
	clk.tic();
	fft.fwd(v);
	double cycles = clk.toc();
	

	MKL_free(space);
	return cycles/count;
}

void printinfo(int m, int n){
	int T[]={1,2,4,8,10,12};
	cout<<setw(50)<<"2D FFT by composing 1D FFTs"<<endl;
	cout<<setw(50)<<"m: "<<m<<endl;
	cout<<setw(50)<<"n: "<<n<<endl;
	cout<<setw(50)<<"nmlzd using division by nm(log2(n)+log2(m))"<<endl;
	double nmlz = m*n*(log10(m*1.0)+log10(n*1.0))/log10(2.0);
	for(int i=5; i < 6; i++){
		double cycles = time2d(m,n,T[i]);
		cout<<endl<<fixed;
		cout<<setw(25)<<"number of threads: "<<T[i]<<endl;
		cout<<setw(25)<<"cycles per FFT (nmlzd): "<<cycles/nmlz<<endl;
		cout<<setw(25)<<"above * T: "<<cycles/nmlz*T[i]<<endl;
	}
	cout<<endl<<endl<<endl;	
}

int main(){
	printinfo(1024, 1024);
	printinfo(27*81, 27*81);
	printinfo(625*5, 625*5);
	printinfo(1024*12, 1024*12);
}
