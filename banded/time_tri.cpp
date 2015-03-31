#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <unistd.h>
#include "../utils/Random.hh"
#include "../utils/utils.hh"
#include "../utils/TimeStamp.hh"
#include "tridiag.hh"
using namespace std;

/*
 * n = matrix dimension
 * nrhs = num of rhs
 * flag = 0 (reuses same solver) (else create new solver every time)
 */
double time_tri(int n, int nrhs, int flag){
	TriSolve tri(n);
	Random rng;
	double *tmp = new double[n];
	for(int i=0; i < n-1; i++)
		tmp[i] = rng.randn();
	tri.setu1(tmp);
	for(int i=0; i < n; i++)
		tmp[i] = rng.randn();
	tri.setdiag(tmp);
	for(int i=0; i < n-1; i++)
		tmp[i] = rng.randn();
	tri.setl1(tmp);

	tri.factor();
	
	int count = 2l*1000*1000*1000/(8*n);
	if(nrhs < 0)
		nrhs = count;
	count = count/nrhs * nrhs;
	double *space = (double *)MKL_malloc(sizeof(double)*n*count, 32);
	double *v = space;
	for(int i=0; i < n*count; i++)
		v[i] = i*i/(i*i+7.0)+1/(i+1.0);
	TimeStamp clk;
	clk.tic();
	for(int i=0; i < count; i+=nrhs){
		if(flag==0)
			tri.solve(v, nrhs);
		else{
			TriSolve trix(n);
			for(int j=0; j < n-1; j++)
				tmp[j] = 1.0/(i+7);
			trix.setu1(tmp);
			for(int j=0; j < n; j++)
				tmp[j] = 10.0/(i+1);
			trix.setdiag(tmp);
			for(int j=0; j < n-1; j++)
				tmp[j] = 1.0/(i+8);
			trix.setl1(tmp);
			trix.factor();
			trix.solve(v, nrhs);
		}
		v += n*nrhs;
	}
	double cycles = clk.toc();
	MKL_free(space);
	delete[] tmp;
	return cycles/count;
}

void print_info(int n, int nrhs, int flag=0){
	cout<<setw(50)<<"Timing Tridiagonal Solver"<<endl;
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
	double cycles = time_tri(n, nrhs, flag);
	cout<<setw(25)<<"cycles: "<<cycles/nmlz<<endl<<endl;
}

int main(){
	int nrhs[] = {1, 2, 10, -1};
	for(int i=0; i < 4; i++){
		print_info(32, nrhs[i]);
		print_info(64, nrhs[i]);
		print_info(128, nrhs[i]);
		print_info(256, nrhs[i]);
		print_info(1024, nrhs[i]);
		print_info(1024*16, nrhs[i]);
		print_info(1024*64, nrhs[i]);
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
