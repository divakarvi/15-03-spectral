#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include "../utils/utils.hh"
#include "../utils/TimeStamp.hh"
#include "../utils/Random.hh"
#include "legendre.hh"
using namespace std;

void cheb(double *restrict x, int M){
	for(int j=0; j <= M; j++)
		x[j] = cos(j*PI/M);
}

//n = size of transform
//fcyc = forward cycles 
//bcyc = backward cycles 
void time_lgndr(int n, double& fcyc, double& bcyc){
	int count = 2l*1000*1000*1000/(8l*(n+1));
	double *space = (double *)MKL_malloc(sizeof(double)*count*(n+1), 32);
	double *v=space;
	FLTrans flt(n);
	TimeStamp clk;
	
	double xx[n+1];
	cheb(xx, n);
	Random rng;

	for(int i=0; i < (n+1)*count; i++){
		double x = xx[i%(n+1)];
		v[i] = pow(x, 3) - x*x;
	}

	
	clk.tic();
	for(int i=0; i < count; i++){
		flt.fwd(v+i*(n+1));
	}
	fcyc = clk.toc()/count;

	clk.tic();
	for(int i=0; i < count; i++){
		flt.bwd(v+i*(n+1));
	}
	bcyc = clk.toc()/count;
	MKL_free(space);
}

void printinfo(int n){
	double fcyc, bcyc;
	time_lgndr(n, fcyc, bcyc);
	cout<<setw(50)<<"Fast Legendre Timing (normalized by n)"<<endl;
	cout<<setw(50)<<"n: "<<n<<endl;
	double nmlz = n*1.0;
	cout<<fixed<<setw(25)<<"forward cycles (nmlz): "<<fcyc/nmlz<<endl;
	cout<<setw(25)<<"backward cycles (nmlz): "<<bcyc/nmlz<<endl;
}

int main(){
	printinfo(64);
	printinfo(128);
	printinfo(256);
	printinfo(512);
	printinfo(1024);
	printinfo(1024*16);
	printinfo(1024*64);
}
