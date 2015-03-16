#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <unistd.h>
#include "../utils/Random.hh"
#include "../utils/utils.hh"
#include "tridiag.hh"
using namespace std;

/*
 * n = matrix dimension
 * leaves data in DBG
 */
void test(int n, int nrhs){
	cout<<setw(50)<<"testing tridiag solve "<<endl;
	cout<<setw(50)<<"n: "<<n<<endl;
	assrt(n>=3);
	char cmd[200];
	sprintf(cmd, "mkdir %s/DBG", getenv("PWD"));
	cout<<setw(10)<<"cmd: "<<cmd<<endl;
	//system(cmd);

	TriSolve tri(n);
	double *tmp = new double[n];
	Random rng;
	for(int i=0; i < n-1; i++)
		tmp[i] = rng.randn();
	tri.setu1(tmp);
	for(int i=0; i < n; i++)
		tmp[i] = rng.randn();
	tri.setdiag(tmp);
	for(int i=0; i < n-1; i++)
		tmp[i] = rng.randn();
	tri.setl1(tmp);
	
	ofstream ofile("DBG/triA.txt");
	ofile<<scientific<<setprecision(16)<<endl;
	double *u1 = tri.getu1();
	double *d = tri.getdiag();
	double *l1 = tri.getl1();
	for(int i=0; i < n; i++)
		if (i==0)
			ofile<<0<<" ";
		else
			ofile<<u1[i-1]<<" ";
	ofile<<endl;
	for(int i=0; i < n; i++)
		ofile<<d[i]<<" ";
	ofile<<endl;
	for(int i=0; i < n; i++)
		if(i<n-1)
			ofile<<l1[i]<<" ";
		else
			ofile<<0<<endl;
	ofile.close();
				

	tri.factor();

	double *f = new double[n*nrhs];
	for(int i=0; i < n*nrhs; i++)
		f[i] = rng.randn();
	array_out(f, n, nrhs, "DBG/trib.txt");
	tri.solve(f, nrhs);
	array_out(f, n, nrhs, "DBG/trix.txt");
	
	cout<<setw(25)<<"test_tridiag.py: "<<endl;
	system("banded.py");
	delete[] tmp;
	delete[] f;
}

int main(){
	test(4, 10);
}
