#include <iostream>
#include <iomanip>
#include <cmath>
#include "../utils/Random.hh"
#include "../utils/utils.hh"
#include "../globals/globals.hh"
#include "../fft/trig.hh"
#include "../pyplot/pyplot.hh"
#include "../diff/diff.hh"
#include "bvp.hh"
#include <omp.h>
using namespace std;

void init_global(int T, int Mi, int Mo){
	gl_init_nthreads(T);
	double pe = 1.53;
	gl_init_grid(Mi, Mo, 100, 100, pe);
	gl_init_boxsize(2, 2);
	gl_init_Re(100);
}

/*
 * (D*D-a*a)u = f with f = -(pi*pi+a*a)*sin(pi*y)
 */
void test1(double a, int Mi, int Mo){
	int T = 1;
	init_global(T, Mi, Mo);
	assrt(gl_state==GL_INIT);
	assrt(gl_Mi == Mi);
	assrt(gl_Mo == Mo);

	int M = Mo*(Mi+1)-1;
	double f[M+1];
	double u[M+1];
	double du[M+1];
	double soln[M+1];
	double dsoln[M+1];
	for(int i=0; i <= M; i++){
		double y = gl_y[i];
		f[i] = -(PI*PI+a*a)*sin(PI*y);
		soln[i] = sin(PI*y);
		dsoln[i] = PI*cos(PI*y);
	}

	BVPSolve bvp(a, Mi, Mo);
	bvp.init_oy(gl_oy);
	bvp.initcheb(gl_chebi);
	DCT dct(Mi);
	bvp.initdct(dct);
	bvp.initbanded();
	bvp.solve(f, u, du, 0.0, 0.0);
	
	{
		//PyPlot plt("u");
		//plt.plot(gl_y, u, M+1);
		//plt.plot(gl_y, soln, M+1);
		//plt.linestyle("--");
		//plt.show();
	}

	array_diff(u, soln, M+1);

	{
		//PyPlot plt("du");
		//plt.plot(gl_y, du, M+1);
		//plt.plot(gl_y, dsoln, M+1);
		//plt.linestyle("--");
		//plt.show();
	}
	array_diff(du, dsoln, M+1);

	cout<<scientific<<setprecision(5)<<endl;
	cout<<"\t\t"<<setw(10)<<"u error: "
	    <<array_max(u, M+1)/array_max(soln, M+1)<<endl;
	cout<<"\t\t"<<setw(10)<<"du error: "
	    <<array_max(du, M+1)/array_max(dsoln, M+1)<<endl;

	gl_exit();
}

void run_test1(){
	double a = 10.0*sqrt(2.0);
	int Mi = 5;
	int Mo = 1000;
	test1(a, Mi, Mo);
}

/*
 * (D*D-a*a)u = f + dg/dy with f = -0.5*(pi*pi+a*a)*sin(pi*y)
 * and g = 0.5*(pi*pi+a*a)*cos(pi*y)/pi
 */
void test1x(double a, int Mi, int Mo){
	int T = 1;
	init_global(T, Mi, Mo);
	assrt(gl_state==GL_INIT);
	assrt(gl_Mi == Mi);
	assrt(gl_Mo == Mo);

	int M = Mo*(Mi+1)-1;
	double f[M+1];
	double g[M+1];
	double u[M+1];
	double du[M+1];
	double soln[M+1];
	double dsoln[M+1];
	for(int i=0; i <= M; i++){
		double y = gl_y[i];
		f[i] = -0.5*(PI*PI+a*a)*sin(PI*y);
		g[i] = 0.5*(PI*PI+a*a)*cos(PI*y)/PI;
		soln[i] = sin(PI*y);
		dsoln[i] = PI*cos(PI*y);
	}

	BVPSolve bvp(a, Mi, Mo);
	bvp.init_oy(gl_oy);
	bvp.initcheb(gl_chebi);
	DCT dct(Mi);
	bvp.initdct(dct);
	bvp.initbanded();
	bvp.solvex(f, g, u, du);
	
	{/*
		PyPlot plt("u");
		plt.plot(gl_y, u, M+1);
		plt.plot(gl_y, soln, M+1);
		plt.linestyle("--");
		plt.show();
	 */}

	array_diff(u, soln, M+1);

	{/*
		PyPlot plt("du");
		plt.plot(gl_y, du, M+1);
		plt.plot(gl_y, dsoln, M+1);
		plt.linestyle("--");
		plt.show();
	 */}
	array_diff(du, dsoln, M+1);

	cout<<scientific<<setprecision(5)<<endl;
	cout<<"\t\t"<<setw(10)<<"u error: "
	    <<array_max(u, M+1)/array_max(soln, M+1)<<endl;
	cout<<"\t\t"<<setw(10)<<"du error: "
	    <<array_max(du, M+1)/array_max(dsoln, M+1)<<endl;

	gl_exit();
}

void run_test1x(){
	double a = 10.0*sqrt(2.0);
	int Mi = 32; //5
	int Mo = 1;  //1000
	test1x(a, Mi, Mo);
}

void init_exact1(double a,
		 double *f, double* soln, double *dsoln, double *ddsoln,
		 double *ylist, int N){
	for(int i=0; i < N; i++){
		double y = ylist[i];
		f[i] = -(PI*PI+a*a)*sin(PI*y);
		soln[i] = sin(PI*y);
		dsoln[i] = PI*cos(PI*y);
		ddsoln[i] = -PI*PI*sin(PI*y);
	}
}

void init_exact2(double a,
		 double *f, double* soln, double *dsoln, double *ddsoln,
		 double *ylist, int N){
	for(int i=0; i < N; i++){
		double y = ylist[i];
		f[i] = 1.0;
		double q1 = exp(-a*(1+y))/(1+exp(-2*a));
		double q2 = exp(a*(-1+y))/(1+exp(-2*a));
		soln[i] = (q1+q2-1)/(a*a);
		dsoln[i] = (-a*q1+a*q2)/(a*a);
		ddsoln[i] = (q1+q2);
	}
}



/*
 * solve (D*D-a*a)u = f 
 * erroru = rel error in u, inf norm
 * errordu1 = rel error in du computed by bvp
 * errordu2 = rel error in du computed by diffing u
 * errorddu1 = rel error in ddu computed by diffing du
 * errorddu2 = rel error in ddu computed by diffing u twice
 */
void errors_bvp(double a, int Mi, int Mo,
		   double& erroru, 
		   double& errordu1, double& errordu2,
		   double& errorddu1, double& errorddu2){
	int T = 1;
	init_global(T, Mi, Mo);
	assrt(gl_state==GL_INIT);
	assrt(gl_Mi == Mi);
	assrt(gl_Mo == Mo);
	
	int N = Mo*(Mi+1);
	double f[N];
	double u[N];
	double du1[N];
	double du2[N];
	double ddu1[N];
	double ddu2[N];
	double soln[N];
	double dsoln[N];
	double ddsoln[N];
	
	init_exact2(a, f, soln, dsoln, ddsoln, gl_y, N);
	
	BVPSolve bvp(a, Mi, Mo);
	bvp.init_oy(gl_oy);
	bvp.initcheb(gl_chebi);
	DCT dct(Mi);
	bvp.initdct(dct);
	bvp.initbanded();
	bvp.solve(f, u, du1, 0.0, 0.0);

	Diff diff(Mi, Mo);
	diff.initcheb(gl_chebi);
	diff.initdct(dct);
	DST dst(Mi);
	diff.initdst(dst);
	diff.init_oy(gl_oy);
	array_copy(du2, u, N);
	diff.diff(du2);
	array_copy(ddu1, du1, N);
	diff.diff(ddu1);
	array_copy(ddu2, du2, N);
	diff.diff(ddu2);

	array_diff(u, soln, N);
	array_diff(du1, dsoln, N);
	array_diff(du2, dsoln, N);
	array_diff(ddu1, ddsoln, N);
	array_diff(ddu2, ddsoln, N);
	
	erroru = array_max(u, N)/array_max(soln, N);
	errordu1 = array_max(du1, N)/array_max(dsoln, N);
	errordu2 = array_max(du2, N)/array_max(dsoln, N);
	errorddu1 = array_max(ddu1, N)/array_max(ddsoln, N);
	errorddu2 = array_max(ddu2, N)/array_max(ddsoln, N);

	cout<<erroru<<endl;
	cout<<errordu1<<endl;
	cout<<errordu2<<endl;
	cout<<errorddu1<<endl;
	cout<<errorddu2<<endl;

	gl_exit();
}

/*
 * test of kleiser-shumann functionality 
 */
void testks(double a, int Mi, int Mo){
	int T = 1;
	init_global(T, Mi, Mo);
	assrt(gl_state==GL_INIT);
	assrt(gl_Mi == Mi);
	assrt(gl_Mo == Mo);

	int M = Mo*(Mi+1)-1;
	double pa[M+1];
	double dpa[M+1];
	double pb[M+1];
	double dpb[M+1];

	BVPSolve bvp(a, Mi, Mo);
	bvp.init_oy(gl_oy);
	bvp.initcheb(gl_chebi);
	DCT dct(Mi);
	bvp.initdct(dct);
	bvp.initbanded();

	bvp.hmg10n01(pa, dpa, pb, dpb);

	double pa_exact[M+1], pb_exact[M+1];
	for(int i=0; i <= M; i++){
		double y = gl_y[i];
		pa_exact[i] = exp(+a*(y-1))/(1-exp(-4*a))
			- exp(-a*(y+1))*exp(-2*a)/(1-exp(-4*a));
		pb_exact[i] =  exp(-a*(y+1))/(1-exp(-4*a))
			- exp(+a*(y-1))*exp(-2*a)/(1-exp(-4*a));
	}

	{
		PyPlot plt("bvperr");
		array_diff(pa, pa_exact, M+1);
		array_diff(pb, pb_exact, M+1);
		array_abs(pa, M+1);
		array_abs(pb, M+1);
		plt.semilogy(gl_y, pa, M+1);
		plt.linewidth("3");
		plt.semilogy(gl_y, pb, M+1);
		plt.title("ks errors (thick line at +1)");
		plt.show();
	}

	{
		PyPlot plt("bvp");
		plt.plot(gl_y, dpa, M+1);
		plt.linewidth("3");
		plt.plot(gl_y, dpb, M+1);
		plt.title("derv of homg solns");
		plt.show();
	}
	gl_exit();
}

int main(){
	//run_test1();
	//run_test1x();

	//double a = 10;
	//int Mi = 32;
	//int Mo = 16;
	//double e1, e2, e3, e4, e5;
	//errors_bvp(a, Mi, 1000, e1, e2, e3, e4, e5);

	double a = 10;
	double Mi = 32;
	double Mo = 34;
	testks(a, Mi, Mo);
}
