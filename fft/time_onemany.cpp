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
 * n = size of c2c FFT
 * count = number of FFTs
 * v = of size 2*n*count
 * fcycles = number of cycles for forward transform
 * bcycles = number of cycles for backard transform
 */
void onecall(int n, int count, double *v,
	     double &fcycles, double& bcycles){
	FFTc2c fft(n, count, n, 1);
	TimeStamp clk;
	clk.tic();
	fft.fwd(v);
	fcycles = clk.toc();
	
	clk.tic();
	fft.bwd(v);
	bcycles = clk.toc();
}

/*
 * n = size of c2c FFT
 * count = number of transforms
 * v = of size 2*n*count
 * fstats = forward fft stats
 * bstats = backward fft stats
 */
void manycalls(int n, int count, double *v, 
	       StatVector& fstats, StatVector& bstats){
	assrt(fstats.getSize() >= count);
	assrt(bstats.getSize() >= count);
	fstats.flush();
	bstats.flush();
	FFTc2c fft(n, 1, n, 1); 
	TimeStamp clk;
	double *p = v;
	double nmlz = n*log(n*1.0)/log(2.0);
	for(int c=0; c < count; c++){
		clk.tic();
		fft.fwd(p);
		fstats.insert(clk.toc()/nmlz);
		p += 2*n;
	}

	for(int c=0; c < count; c++){
		p -= 2*n;
		clk.tic();
		fft.bwd(p);
		bstats.insert(clk.toc()/nmlz);
	}
}

void printinfo(int n){
	cout<<setw(50)<<"******************************"<<endl;
	cout<<setw(50)<<"ALL OUTPUT cycles/nlg2n "<<endl;
	cout<<setw(50)<<"******************************"<<endl<<endl;
	int count = 8l*1000*1000*1000/(sizeof(double)*4l*n); //8 GB of space
	double *space = (double *)MKL_malloc(sizeof(double)*4*n*count, 32);
	double *v = space;
	double *w = v + 2*n*count;
	for(long i=0; i < 2*n*count; i++)
		v[i] = w[i] = 1.0/(i+1);

	cout<<setw(50)<<"n for c2c fft: "<<n<<endl;
	cout<<setw(50)<<"number of transforms: "<<count<<endl;
	cout<<endl<<setw(50)<<"single call many transforms"<<endl;
	double fcycles;
	double bcycles;
	onecall(n, count, v, fcycles, bcycles);
	double nmlz = 1.0*n*log(n*1.0)/log(2.0);
	cout<<setw(25)<<"avg cycles per fwd fft:"<<fcycles/count/nmlz<<endl;
	cout<<setw(25)<<"avg cycles per bwd fft:"<<bcycles/count/nmlz<<endl;

	cout<<endl<<setw(50)<<"many calls many transforms"<<endl;
	StatVector fstats(count), bstats(count);
	manycalls(n, count, v, fstats, bstats);
	fstats.print("forward transform");
	bstats.print("backward transform");

	array_diff(v, w, 2*n*count);
	cout<<endl<<setw(25)<<"ERROR: "
	    <<array_max(v, 2*n*count)/array_max(w, 2*n*count)<<endl<<endl;

	MKL_free(space);
}


int main(){
	printinfo(256);
	printinfo(1024);
	printinfo(1024*12);
}
