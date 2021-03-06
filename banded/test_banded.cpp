#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <unistd.h>
#include "../utils/Random.hh"
#include "../utils/utils.hh"
#include "banded.hh"
using namespace std;

/*
 * n = matrix dimension
 * leaves data in DBG
 */
void banded(int n, int nrhs){
	cout<<setw(50)<<"testing banded solve with l=u=2"<<endl;
	cout<<setw(50)<<"n: "<<n<<endl;
	assrt(n>=4);
	char cmd[200];
	sprintf(cmd, "mkdir %s/DBG", getenv("PWD"));
	cout<<setw(10)<<"cmd: "<<cmd<<endl;
	//system(cmd);

	BandedSolve banded(n);
	double *tmp = new double[n*nrhs];
	Random rng;
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
	double *A = banded.getA();
	int lda = 7;
	array_out(A+2, 5, n, lda, "DBG/banded_A.txt");//must precede factor()
 
	banded.factor();

	for(int i=0; i < n*nrhs; i++)
		tmp[i] = rng.randn();
	array_out(tmp, n, nrhs, "DBG/b.txt");
	banded.solve(tmp, nrhs);
	array_out(tmp, n, nrhs, "DBG/x.txt");

	cout<<setw(25)<<"calling banded.py: "<<endl;
	system("banded.py 2 2");

	delete[] tmp;
}

void cholesky(){
	int n = 5;
	double u2[n-2], u1[n-1], d[n];
	for(int j=0; j < n-2; j++)
		u2[j] = 0;
	for(int j=0; j < n-1; j++)
		u1[j] = -1;
	for(int j=0; j < n; j++)
		d[j] = 2.0;

	CholeskyBanded cb(n);
	cb.setu2(u2);
	cb.setu1(u1);
	cb.setdiag(d);
	cb.factor();

	double *A = cb.getA();
	array_out(A, 3, n, 3, "tmp.dat");
}

int main(){
	//banded(5, 10);
	cholesky();
}
