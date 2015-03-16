#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include "../utils/Random.hh"
#include "../utils/utils.hh"
#include "trig.hh"
using namespace std;

/*
 * test Discrete Cosine Transform with simple harmonics
 */
void test_dct(int n){
	double *space;
	space = (double *)MKL_malloc(sizeof(double)*(n+1)*2, 32);
	double *v=space;
	double *w = v + (n+1);
	assrt(n>=7);
	for(int i=0; i <= n; i++){
		double x = 1.0*i*PI/n;
		v[i] = 2.0 - 11.0*cos(x) + 7.0*cos(7*x);
		w[i] = v[i];
	}

	cout<<setw(50)<<"DCT of 2-11cos(x)+7cos(7x)"<<endl;
	cout<<setw(50)<<"n: "<<n<<endl;
	cout<<fixed;
	DCT dct(n);
	dct.fwd(v);
	for(int i=0; i <= 7; i++)
		cout<<setw(25)<<"coeff = "<<v[i]<<endl;
	dct.bwd(v);
	array_diff(v, w, n+1);
	double error = array_max(v, n+1)/array_max(w, n+1);
	cout<<setw(25)<<"Error: "<<scientific<<error<<endl<<endl;
	MKL_free(space);
}

void test_dst(int n){
	double *space;
	space = (double *)MKL_malloc(sizeof(double)*(n+1)*2, 32);
	double *v=space;
	double *w = v + (n+1);
	assrt(n>=7);
	for(int i=0; i <= n; i++){
		double x = 1.0*i*PI/n;
		v[i] = - 11.0*sin(x) + 7.0*sin(7*x);
		w[i] = v[i];
	}

	cout<<setw(50)<<"DST of -11sin(x)+7sin(7x)"<<endl;
	cout<<setw(50)<<"n: "<<n<<endl;
	cout<<fixed;
	DST dst(n);
	dst.fwd(v);
	for(int i=0; i <= 7; i++)
		cout<<setw(25)<<"coeff = "<<v[i]<<endl;
	dst.bwd(v);
	array_diff(v, w, n+1);
	double error = array_max(v, n+1)/array_max(w, n+1);
	cout<<setw(25)<<"Error: "<<scientific<<error<<endl<<endl;
	MKL_free(space);
}

int main(){
	test_dct(7);
	test_dct(8);
	test_dct(128);

	test_dst(7);
	test_dst(8);
	test_dst(128);
}
