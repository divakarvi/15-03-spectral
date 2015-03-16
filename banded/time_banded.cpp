#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <unistd.h>
#include "../utils/Random.hh"
#include "../utils/utils.hh"
#include "../utils/TimeStamp.hh"
#include "banded.hh"
using namespace std;

/*
 * n = matrix dimension
 * nrhs = num of rhs
 * flag = 0 (reuses same solver) (else create new solver every time)
 */
double time_banded(int n, int nrhs, int flag){
	BandedSolve banded(n);
	Random rng;
	double *tmp = new double[n];
	for(int i=0; i < n-2; i++)
		tmp[i] = rng.randn();
	banded.setu2(tmp);
	for(int i=0; i < n-1; i++)
		tmp[i] = rng.randn();
	banded.setu1(tmp);
	for(int i=0; i < n; i++)
		tmp[i] = rng.randn();
	banded.setdiag(tmp);
	for(int i=0; i < n-1; i++)
		tmp[i] = rng.randn();
	banded.setl1(tmp);
	for(int i=0; i < n-2; i++)
		tmp[i] = rng.randn();
	banded.setl2(tmp);

	banded.factor();
	
	int count = 4l*1000*1000*1000/(8*n);
	if(nrhs < 0)
		nrhs = count;
	count = count/nrhs * nrhs;
	double *space = (double *)MKL_malloc(sizeof(double)*n*count, 32);
	double *v = space;
	for(int i=0; i < n*count; i++)
		v[i] = i*i/(i*i+6.0)+1/(i+2.0);
	TimeStamp clk;
	clk.tic();
	for(int i=0; i < count; i+=nrhs){
		if(flag==0)
			banded.solve(v, nrhs);
		else{
			BandedSolve bandedx(n);
			for(int j=0; j < n-1; j++)
				tmp[j] = 1.0/(i+20);
			bandedx.setu2(tmp);
			bandedx.setu1(tmp);
			for(int j=0; j < n; j++)
				tmp[j] = 20.0/(i+1);
			bandedx.setdiag(tmp);
			for(int j=0; j < n-1; j++)
				tmp[j] = 1.0/(i+8);
			bandedx.setl1(tmp);
			bandedx.setl2(tmp);
			bandedx.factor();
			bandedx.solve(v, nrhs);
		}
		v += n*nrhs;
	}
	double cycles = clk.toc();
	MKL_free(space);
	delete[] tmp;
	return cycles/count;
}

void print_info(int n, int nrhs, int flag=0){
	cout<<setw(50)<<"Timing Banded (l,u =2) Solver"<<endl;
	if(flag==0)
		cout<<setw(50)<<"reuses same factrization"<<endl;
	else
		cout<<setw(50)<<"new matrix every trial"<<endl;
	cout<<setw(50)<<"n: "<<n<<endl;
	if(nrhs > 0)
		cout<<setw(50)<<"nrhs: "<<nrhs<<endl;
	else
		cout<<setw(50)<<"nrhs: max possible"<<endl;
	cout<<setw(50)<<"nmlzed by n"<<endl;
	double nmlz = n;
	double cycles = time_banded(n, nrhs, flag);
	cout<<setw(25)<<"cycles: "<<cycles/nmlz<<endl<<endl;
}


int main(){
	int nrhs[] = {1, 2, 10, -1};
	for(int i=0; i < 4; i++){
		print_info(16, nrhs[i]);
		print_info(32, nrhs[i]);
		print_info(64, nrhs[i]);
		print_info(128, nrhs[i]);
		print_info(256, nrhs[i]);
	}
	

	cout<<endl<<setw(50)<<"*** new matrix every trial ***"<<endl<<endl;
	
	for(int i=0; i < 4; i++){
		print_info(16, nrhs[i], -1);
		print_info(32, nrhs[i], -1);
		print_info(64, nrhs[i], -1);
		print_info(128, nrhs[i], -1);
		print_info(256, nrhs[i], -1);
	}
}
	

			
