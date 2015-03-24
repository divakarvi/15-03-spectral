#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <unistd.h>
#include <omp.h>
#include "../utils/Random.hh"
#include "../utils/utils.hh"
#include "../banded/tridiag.hh"
#include "../fft/trig.hh"
#include "bvpi.hh"
using namespace std;

void test_hmg(double a, int M){
	BVPSolvei bvpi(a, M);
	double chebpts[M+1];
	cheb(chebpts, M);
	bvpi.initcheb(chebpts);
	DCT dct(M);
	bvpi.initdct(dct);
	bvpi.inithmg();

	cout<<setw(50)<<"testing homg solns of (D^2-a^2)u = 0"<<endl;
	cout<<setw(50)<<"a: "<<a<<endl;
	cout<<setw(50)<<"M: "<<M<<endl<<endl;

	double *h1 = bvpi.geth1();
	double *d1 = bvpi.getd1();
	dct.fwd(h1);
	dct.fwd(d1);
	array_out(&a, 1, 1, "DBG/a.txt");
	array_out(h1, M+1, 1, "DBG/u.txt");
	array_out(d1, M+1, 1, "DBG/du.txt");
	for(int i=0; i <= M; i++)
		h1[i] = 0;
	array_out(h1, M+1, 1, "DBG/f.txt");

	cout<<setw(50)<<"verifying T_0 u = 1, T_0 Du = 0"<<endl;
	system("bvpi_verify.py");
	cout<<endl;

	double *h2 = bvpi.geth2();
	double *d2 = bvpi.getd2();
	dct.fwd(h2);
	dct.fwd(d2);
	array_out(h2, M+1, 1, "DBG/u.txt");
	array_out(d2, M+1, 1, "DBG/du.txt");
	
	cout<<setw(50)<<"verifying T_0 u = 0, T_0 Du = 1"<<endl;
	system("bvpi_verify.py");
	cout<<endl;
}


/*
 * aa[0..T-1] = a value to be used by thread
 * M = dimension of bvpi
 * nrhs = nrhs per thread
 * T = number of threads 
 * test with g = 0
 */
void test_bvpi(double *aa, int M, int nrhs, int T){
	double gl_chebi[M+1];
	cheb(gl_chebi, M);
	double *f = (double *)MKL_malloc((M+1)*nrhs*T*sizeof(double), 64);
	double *u = (double *)MKL_malloc((M+1)*nrhs*T*sizeof(double), 64);
	double *du = (double *)MKL_malloc((M+1)*nrhs*T*sizeof(double), 64);
	
#pragma omp parallel				\
	num_threads(T)				\
	default(none)				\
	shared(aa, M, nrhs, f, u, du, gl_chebi)
	{
		int tid = omp_get_thread_num();
		double a = aa[tid];
		DCT dct(M);
		BVPSolvei bvpi(a, M);
		bvpi.initcheb(gl_chebi);
		bvpi.initdct(dct);
		bvpi.inithmg();
		double *h1 = bvpi.geth1();
		double *h2 = bvpi.geth2();
		double *d1 = bvpi.getd1();
		double *d2 = bvpi.getd2();

		double *ff = f + tid*(M+1)*nrhs;
		double *uu = u + tid*(M+1)*nrhs;
		double *ddu = du + tid*(M+1)*nrhs;
		for(int i=0; i < nrhs; i++)
			for(int j=0; j <= M; j++){
				double y = gl_chebi[j];
				ff[j+i*(M+1)] = -(PI*PI+a*a)*sin(PI*y);
			}
		bvpi.solvep(ff, uu, ddu, nrhs);
		double B[4];
		B[0] = -h1[0];
		B[1] = h1[M];
		B[2] = -h2[0];
		B[3] = h2[M];
		double rhs[2];
		double e[2];
		for(int i=0; i < nrhs; i++){
			rhs[0] = uu[i*(M+1)];
			rhs[1] = -uu[M+i*(M+1)];
			solve2x2(B, rhs, e);
			for(int j=0; j <= M; j++){
				uu[j+i*(M+1)] += e[0]*h1[j] + e[1]*h2[j];
				ddu[j+i*(M+1)] += e[0]*d1[j] + e[1]*d2[j];
			}
		}
	}//end parallel region
	
	DCT dct(M);

	cout<<"\t\t\t\t"<<"Verifying BVPSolvei"<<endl;
	cout<<"\t\t\t\t"<<"(D^2-a^2)u = -(pi*pi+a*a)*sin(pi*y)"<<endl<<endl;
	for(int t=0; t < T; t++)
		for(int i=0; i < nrhs; i++){
			double a = aa[t];
			double *ff = f + i*(M+1)+t*(M+1)*nrhs;
			double *uu = u + i*(M+1)+t*(M+1)*nrhs;
			double *ddu = du + i*(M+1)+t*(M+1)*nrhs;
			cout<<endl;
			cout<<"\t\t"<<setw(10)<<"a: "<<a<<endl;
			cout<<"\t\t"<<setw(10)<<"i: "<<i<<endl;
			cout<<"\t\t"<<setw(10)<<"u at 1: "<<uu[0]<<endl;
			cout<<"\t\t"<<setw(10)<<"u at -1: "<<uu[M]<<endl;
			dct.fwd(ff);
			dct.fwd(uu);
			dct.fwd(ddu);
			array_out(&a, 1, 1, "DBG/a.txt");
			array_out(uu, M+1, 1, "DBG/u.txt");
			array_out(ddu, M+1, 1, "DBG/du.txt");
			array_out(ff, M+1, 1, "DBG/f.txt");
			/*
			 * bndries not verified---I think
			 */
			system("bvpi_verify.py");
		}

	MKL_free(f);
	MKL_free(u);
	MKL_free(du);
}


/*
 * aa[0..T-1] = a value to be used by thread
 * M = dimension of bvpi
 * nrhs = nrhs per thread
 * T = number of threads 
 * test with g != 0
 */
void test_bvpix(double *aa, int M, int nrhs, int T){
	double gl_chebi[M+1];
	cheb(gl_chebi, M);
	double *f = (double *)MKL_malloc((M+1)*nrhs*T*sizeof(double), 64);
	double *g = (double *)MKL_malloc((M+1)*nrhs*T*sizeof(double), 64);
	double *u = (double *)MKL_malloc((M+1)*nrhs*T*sizeof(double), 64);
	double *du = (double *)MKL_malloc((M+1)*nrhs*T*sizeof(double), 64);
	
#pragma omp parallel				\
	num_threads(T)				\
	default(none)				\
	shared(aa, M, nrhs, f, g, u, du, gl_chebi)
	{
		int tid = omp_get_thread_num();
		double a = aa[tid];
		DCT dct(M);
		BVPSolvei bvpi(a, M);
		bvpi.initcheb(gl_chebi);
		bvpi.initdct(dct);
		bvpi.inithmg();
		double *h1 = bvpi.geth1();
		double *h2 = bvpi.geth2();
		double *d1 = bvpi.getd1();
		double *d2 = bvpi.getd2();

		double *ff = f + tid*(M+1)*nrhs;
		double *gg = g + tid*(M+1)*nrhs;
		double *uu = u + tid*(M+1)*nrhs;
		double *ddu = du + tid*(M+1)*nrhs;
		for(int i=0; i < nrhs; i++)
			for(int j=0; j <= M; j++){
				double y = gl_chebi[j];
				ff[j+i*(M+1)] = -0.5*(PI*PI+a*a)*sin(PI*y);
				gg[j+i*(M+1)] = 0.5*(PI*PI+a*a)*cos(PI*y)/PI;
			}
		bvpi.solvepx(ff, gg, uu, ddu, nrhs);
		double B[4];
		B[0] = -h1[0];
		B[1] = h1[M];
		B[2] = -h2[0];
		B[3] = h2[M];
		double rhs[2];
		double e[2];
		for(int i=0; i < nrhs; i++){
			rhs[0] = uu[i*(M+1)];
			rhs[1] = -uu[M+i*(M+1)];
			solve2x2(B, rhs, e);
			for(int j=0; j <= M; j++){
				uu[j+i*(M+1)] += e[0]*h1[j] + e[1]*h2[j];
				ddu[j+i*(M+1)] += e[0]*d1[j] + e[1]*d2[j];
			}
		}
	}//end parallel region
	
	DCT dct(M);

	cout<<"\t\t\t\t"<<"Verifying BVPSolvei"<<endl;
	cout<<"\t\t\t\t"<<"(D^2-a^2)u = -(pi*pi+a*a)*sin(pi*y)"<<endl<<endl;
	for(int t=0; t < T; t++)
		for(int i=0; i < nrhs; i++){
			double a = aa[t];
			double *ff = f + i*(M+1)+t*(M+1)*nrhs;
			for(int j=0; j <= M; j++)
				ff[j] *= 2.0;
			double *uu = u + i*(M+1)+t*(M+1)*nrhs;
			double *ddu = du + i*(M+1)+t*(M+1)*nrhs;
			cout<<endl;
			cout<<"\t\t"<<setw(10)<<"a: "<<a<<endl;
			cout<<"\t\t"<<setw(10)<<"i: "<<i<<endl;
			cout<<"\t\t"<<setw(10)<<"u at 1: "<<uu[0]<<endl;
			cout<<"\t\t"<<setw(10)<<"u at -1: "<<uu[M]<<endl;
			dct.fwd(ff);
			dct.fwd(uu);
			dct.fwd(ddu);
			array_out(&a, 1, 1, "DBG/a.txt");
			array_out(uu, M+1, 1, "DBG/u.txt");
			array_out(ddu, M+1, 1, "DBG/du.txt");
			array_out(ff, M+1, 1, "DBG/f.txt");
			/*
			 * bndries not verified---I think
			 */
			system("bvpi_verify.py");
		}

	MKL_free(f);
	MKL_free(u);
	MKL_free(du);
}



int main(){
	verify_dir("DBG");
	{/*
		double a[] = {1, 20, 47.7, .01};
		int M[] = {20, 50, 100, 5};
		for(int i=0; i < 4; i++)
			test_hmg(a[i], M[i]);
	 */}

	{/*
		const int T = 16;
		double aa[T];
		for(int t=0; t < T; t++)
			aa[t] = (t+2)*(t+3)*(t+4)*(t+5)+sqrt(PI); 
		int M = 42;
		int nrhs = 4;
		test_bvpi(aa, M, nrhs, T);
	 */}

	{
		const int T = 16;
		double aa[T];
		for(int t=0; t < T; t++)
			aa[t] = (t+2)*(t+3)*(t+4)*(t+5)+sqrt(PI); 
		int M = 42;
		int nrhs = 4;
		test_bvpix(aa, M, nrhs, T);
	}
}
