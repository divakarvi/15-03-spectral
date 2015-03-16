#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include "../utils/Random.hh"
#include "../utils/utils.hh"
#include "../utils/TimeStamp.hh"
#include "../utils/StatVector.hh"
#include "fft_mkl.hh"
using namespace std;

/*
 * n = size of fft (1D c2c) (unit = complex)
 * stride (unit = complex)
 * return value = avg number of cycles for a single fwd trans of size 
 */
double  time_stride(int n, int stride){
	double *space;
	int sset = 2*8*stride*n; //unit = byte
	int count = 2l*1000*1000*1000/sset;
	space = (double *)MKL_malloc(sset*count, 32);
	double *v = space;
	for(long i=0; i < 2*stride*n*count; i++)
		v[i] = 1.0/(i*i+7);
	
	TimeStamp clk;
	FFTc2c fft(n, stride, 1, stride);
	clk.tic();
	for(int i=0; i < count; i++){
		fft.fwd(v);
		v += 2*n*stride;
	}
	double cycles = clk.toc();

	MKL_free(space);

	return cycles/stride/count;
}

int main(){
	cout<<setw(50)<<"******************************"<<endl;
	cout<<setw(50)<<"ALL OUTPUT cycles/nlg2n "<<endl;
	cout<<setw(50)<<"******************************"<<endl<<endl;
	int n = 1024;
	double cyc;
	for(int k=1; k < 9; k+=7)
		for(int stride=k*n; stride < k*n+4; stride++){
			cyc = time_stride(n, stride);
			double nmlz = n*log(1.0*n)/log(2.0);
			cout<<setw(25)<<"n: "<<n<<endl;
			cout<<setw(25)<<"stride: "<<stride<<endl;
			cout<<setw(25)<<"cycles: "<<cyc/nmlz<<endl<<endl;
		}
}
