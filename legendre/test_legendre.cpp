#include <cmath>
#include "../utils/utils.hh"
#include "../utils/Random.hh"
#include "../pyplot/pyplot.hh"
#include "../fft/trig.hh"
#include "legendre.hh"

void cheb(double *restrict x, int M){
	for(int j=0; j <= M; j++)
		x[j] = cos(j*PI/M);
}

void test1(){
	int M = 64;
	double x[M+1], f[M+1];
	cheb(x, M);
	
	/*
	 * f = 1 + P_3(x)
	 */
	for(int i=0; i <= M; i++)
		f[i] = 1 + 0.5*(5*pow(x[i], 3) - 3*x[i]);
	PyPlot plt("P1P3");
	plt.plot(x, f, M+1);
	plt.title("legendre poly 1 + P3");
		
	FLTrans flt(M);
	flt.fwd(f);
	array_show(f, 10, "legendre coeffs");
}

void test2(){
	int M = 64;
	double x[M+1], f[M+1], ff[M+1];
	cheb(x, M);
	
	/*
	 * f = 1 + 2*P_1 + 4*P_2 + 8*P_3
	 */
	for(int i=0; i <= M; i++)
		ff[i] = f[i] = 1 + 2*x[i] 
			+ 4*(0.5*(3*x[i]*x[i]-1)) 
			+ 8*(0.5*(5*pow(x[i], 3) - 3*x[i]));
		
	FLTrans flt(M);
	flt.fwd(f);
	array_show(f, 10, "legendre coeffs");

	flt.bwd(f);
	array_diff(ff, f, M+1);
	array_abs(ff, M+1);

	PyPlot plt("error");
	plt.semilogy(ff, M+1);
	plt.linestyle("None");
	plt.markertype("o");
	plt.markercolor("black");
	plt.title("fwd/bwd abs error");
	plt.show();
}

void test3(int M){
	double f[M+1], ff[M+1], x[M+1];
	Random rng;
	cheb(x, M);
	for(int i=0; i <= M; i++)
		ff[i] = f[i] = 1.0 + x[i] + x[i]*x[i]+1e-6*rng.randn();

	FLTrans flt(M);
	array_show(f, 10, "f");
	flt.fwd(f);
	flt.bwd(f);
	array_show(f, 10, "f after fwd/bwd");
	
	array_diff(ff, f, M+1);
	array_abs(ff, M+1);

	PyPlot plt("error");
	plt.semilogy(ff, M+1);
	plt.linestyle("None");
	plt.markertype("o");
	plt.markercolor("black");
	plt.title("fwd/bwd error, f = 1+x+x^2+1e-6*noise");
	plt.show();
}

int main(){
	//test1();
	//test2();
	test3(1024);
}
