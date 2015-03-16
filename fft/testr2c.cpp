#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include "../utils/Random.hh"
#include "../utils/utils.hh"
#include "fft_mkl.hh"
using namespace std;

/*
 * single vector
 * superposed harmonics
 * check Fourier coeffs
 */
void testr2c_1(int n){
	double *space;
	space = (double *)MKL_malloc(sizeof(double)*4*(n/2+1), 32);
	double *v = space;
	double *w = v+2*(n/2+1);
	for(int i=0; i < n; i++){
		double x = 2*PI*i/n;
		v[i] = 2.0*cos(x)-4.0*sin(3*x);
		w[i] = v[i];
	}
	
	FFTr2c fft(n, 1, 2*(n/2+1));
	fft.fwd(v);
	cout<<setw(50)<<"r2c FFT"<<endl;
	cout<<setw(50)<<"n: "<<n<<endl;
	cout<<setw(50)<<"fn: "<<"2cos(x)-4sin(3x)"<<endl;
	assrt(n>18);
	for(int i=0; i < 8; i++)
		cout<<fixed<<setw(25)<<"Fourier Coeff:"
		    <<"("<<v[2*i]<<","<<v[2*i+1]<<")"<<endl;
	cout<<endl<<endl;
	cout<<setw(50)<<"Two unused spots in CCE format garbled"<<endl;
	v[1] = 1e100;//garble unused sites
	v[n+1] = -1e100; 
	fft.bwd(v);
	array_diff(v, w, n);
	cout<<scientific<<setw(25)<<"error: "
	    <<array_max(v, n)/array_max(w, n)
	    <<endl;

	MKL_free(space);
}

void testr2c_2(int n, int count, int dist){
	assrt(dist>=2*(n/2+1));
	double *space;
	space = (double *)MKL_malloc(sizeof(double)*2*dist*count, 32);
	double *v = space;
	double *w = v+dist*count;
	for(int i=0; i < dist*count; i++){
		int j = i%dist;
		int blk = i/dist;
		double x = 2.0*PI*j/n;
		if(j>n)
			v[i] = -1;
		else
			v[i] = (blk+1)*2.0*cos(x)-(blk+1)*4.0*sin(3*x);
		w[i] = v[i];
	}
	
	cout<<setw(50)<<"r2c FFT test with superposed harmonics"<<endl;
	cout<<setw(50)<<"n: "<<n<<endl;
	cout<<setw(50)<<"number of blocks: "<<count<<endl;
	cout<<setw(50)<<"distance between blocks: "<<dist<<endl;
	cout<<setw(50)<<"block b: "<<"2(b+1)*cos(x)-4i*(b+1)*sin(3x)"<<endl;
	cout<<setw(25)<<"inf norm of v: "<<array_max(v, dist*count)<<endl;
	FFTr2c fft(n, count, dist);
	fft.fwd(v);
	int blk = count/2;
	cout<<setw(25)<<"block: "<<blk<<endl;
	for(int i=0; i < 8; i++)
		cout<<setw(25)<<"Fourier coeff blk : "
		    <<"("<<v[blk*dist+2*i]<<","<<v[blk*dist+2*i+1]<<")"<<endl;
	fft.bwd(v);
	array_diff(v, w, dist*count);
	double verr=0, winf=0;
	for(int i=0; i < count; i++)
		for(int j=0; j < n; j++){
			verr += fabs(v[i*dist+j]);
			winf += fabs(w[i*dist+j]);
		}
	cout<<endl<<setw(25)<<"n: "<<setw(25)<<n<<endl;
	cout<<setw(25)<<"L1 error: "<<setw(25)<<scientific
	    <<verr/winf<<endl;
	
	MKL_free(space);
}

int main(){
	testr2c_1(2056);
	testr2c_2(1024, 128, 1028); 
}
