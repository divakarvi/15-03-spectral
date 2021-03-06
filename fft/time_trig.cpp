#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include "../utils/Random.hh"
#include "../utils/utils.hh"
#include "../utils/TimeStamp.hh"
#include "trig.hh"
using namespace std;

//n = size of transform
//fcyc = forward cycles 
//bcyc = backward cycles 
void time_dct(int n, double& fcyc, double& bcyc){
	int count = 2l*1000*1000*1000/(8l*(n+1));
	double *space = (double *)MKL_malloc(sizeof(double)*count*(n+1), 32);
	double *v=space;
	DCT dct(n);
	TimeStamp clk;
	
	for(int i=0; i < (n+1)*count; i++)
		v[i] = i*i/(i*i+7.0)+1/(i+1.0);

	clk.tic();
	for(int i=0; i < count; i++){
		dct.fwd(v+i*(n+1));
	}
	fcyc = clk.toc()/count;

	clk.tic();
	for(int i=0; i < count; i++){
		dct.bwd(v+i*(n+1));
	}
	bcyc = clk.toc()/count;
	MKL_free(space);
}

void printinfo(int n){
	double fcyc, bcyc;
	time_dct(n, fcyc, bcyc);
	cout<<setw(50)<<"DCT Timing (normalized by n)"<<endl;
	cout<<setw(50)<<"n: "<<n<<endl;
	double nmlz = n*1.0;
	cout<<fixed<<setw(25)<<"forward cycles (nmlz): "<<fcyc/nmlz<<endl;
	cout<<setw(25)<<"backward cycles (nmlz): "<<bcyc/nmlz<<endl;
}

int main(){
	printinfo(64);
	printinfo(128);
	printinfo(256);
	printinfo(1024);
	printinfo(1024*16);
	printinfo(1024*64);
}
