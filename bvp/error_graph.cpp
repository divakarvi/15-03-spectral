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

static const char *l_title = "no title";
/*
 * instances of (D*D-a*a)u = f
 * SINPIY: u = sin(pi*y)
 * SINKPIY: u = sin(k*pi*y) k adjusted to _6_ points per wavelength (k=Mi*Mo/6)
 * FHYY: u = sin(pi*y)/(1+400y*y) (Four Hundred YY)
 * ONE: f = 1
 */
enum bvp_ex {SINPIY, SINKPIY, FHYY, ONE, NSMTH};

class BVPErrors{
private:
	int state;
	double a;
	int Mi; 
	int Mo;
	BVPSolve bvp;
	DCT dct;
	DST dst;
	Diff diff;
	double *y;
	double *f, *soln, *dsoln, *ddsoln;
	double *u, *du1, *du2, *ddu1, *ddu2;
	double erroru, errordu1, errordu2, errorddu1, errorddu2;
public:
	BVPErrors(double aa, int Minner, int Mouter);
	~BVPErrors();
	void inity(double *yi);
	/*
	 * see above for flag
	 */
	void initexact(enum bvp_ex flag);
	void initcheb(const double *oy, double *chebi);
	void computesolns();
	void computeerrors();
	double error_u(){assrt(state==6); return erroru;}
	double error_du1(){assrt(state==6); return errordu1;}
	double error_du2(){assrt(state==6); return errordu2;}
	double error_ddu1(){assrt(state==6); return errorddu1;}
	double error_ddu2(){assrt(state==6); return errorddu2;}
	double *getsoln();
	double *getdsoln();
	double *getddsoln();
};

BVPErrors::BVPErrors(double aa, int Minner, int Mouter)
	:state(0), a(aa), Mi(Minner), Mo(Mouter),
	 bvp(a, Mi, Mo), dct(Mi), dst(Mi), diff(Mi, Mo)
{
	state = 1;
	int N = (Mi+1)*Mo;
	f = new double[N];
	soln = new double[N];
	dsoln = new double[N];
	ddsoln = new double[N];
	u = new double[N];
	du1 = new double[N];
	du2 = new double[N];
	ddu1 = new double[N];
	ddu2 = new double[N];
}

BVPErrors::~BVPErrors(){
	delete f;
	delete soln;
	delete dsoln;
	delete ddsoln;
	delete u;
	delete du1;
	delete du2;
	delete ddu1;
	delete ddu2;
}

void BVPErrors::inity(double *yin){
	assrt(state==1);
	state = 2;
	y = yin;
}

/*
 * instances of (D*D-a*a)u = f
 * SINPIY: u = sin(pi*y)
 * FHYY: u = sin(pi*y)/(1+400y*y)
 * ONE: f = 1
 * NSMTH: 2nd derv of u has jump at 0.25
 */
void BVPErrors::initexact(enum bvp_ex FLAG){
	assrt(state==2);
	state = 3;
	int N = (Mi+1)*Mo;
	Random rng;
	int k;
	switch(FLAG){
	case SINPIY:
		k = 1;
		for(int i=0; i < N; i++){
			f[i] = -(k*k*PI*PI+a*a)*sin(k*PI*y[i]);
			soln[i] = sin(k*PI*y[i]);
			dsoln[i] = k*PI*cos(k*PI*y[i]);
			ddsoln[i] = -k*k*PI*PI*sin(k*PI*y[i]);
		}	
		break;
	case SINKPIY:
		k = Mi*Mo/6;
		for(int i=0; i < N; i++){
			f[i] = -(k*k*PI*PI+a*a)*sin(k*PI*y[i]);
			soln[i] = sin(k*PI*y[i]);
			dsoln[i] = k*PI*cos(k*PI*y[i]);
			ddsoln[i] = -k*k*PI*PI*sin(k*PI*y[i]);
		}	
		break;
	case FHYY:
		k = 20;
		for(int i=0; i < N; i++){
			double yy  = y[i];
			soln[i] = sin(PI*yy)/(1+k*k*yy*yy);
			dsoln[i] = PI*cos(PI*yy)/(1+k*k*yy*yy)
				- sin(PI*yy)*2*k*k*yy/pow(1+k*k*yy*yy,2.0);
			ddsoln[i] = -PI*PI*sin(PI*yy)/(1+k*k*yy*yy)
				-4*PI*cos(PI*yy)*k*k*yy/pow(1+k*k*yy*yy,2.0)
				+sin(PI*yy)*(6.0*k*k*k*k*yy*yy-2.0*k*k)
				*pow(1+k*k*yy*yy, -3.0);
			f[i] = ddsoln[i]-a*a*soln[i];
		}
		break;
	case ONE:
		for(int i=0; i < N; i++){
			f[i] = 1.0;
			double q1 = exp(-a*(1+y[i]))/(1+exp(-2*a));
			double q2 = exp(a*(-1+y[i]))/(1+exp(-2*a));
			soln[i] = (q1+q2-1)/(a*a);
			dsoln[i] = (-a*q1+a*q2)/(a*a);
			ddsoln[i] = (q1+q2);
		}
		break;
	case NSMTH:
		for(int i=0; i < N; i++){
			int s = (y[i]>0.25)?1:-1;
			double c = (y[i]-0.25);
			double b = (y[i]+1);
			soln[i] = s*c*c*c/6.0+(49.0/384)*b-125.0/384;
			dsoln[i] = s*c*c/2.0 +(49.0/384);
			ddsoln[i] = s*c;
			f[i] = ddsoln[i] - a*a*soln[i];
		}
		break;
	}
}

void BVPErrors::initcheb(const double *oy, double *chebi){
	assrt(state==3);
	state = 4;
	bvp.init_oy(oy);
	bvp.initcheb(chebi);
	bvp.initdct(dct);			
	bvp.initbanded();
	
	diff.initcheb(chebi);
	diff.initdct(dct);
	diff.initdst(dst);
	diff.init_oy(oy);
}

void BVPErrors::computesolns(){
	assrt(state==4);
	state = 5;
	int N = (Mi+1)*Mo;
	bvp.solve(f, u, du1);
	array_copy(du2, u, N);
	diff.diff(du2);

	array_copy(ddu1, du1, N);
	diff.diff(ddu1);
	array_copy(ddu2, du2, N);
	diff.diff(ddu2);
}

void BVPErrors::computeerrors(){
	assrt(state==5);
	state = 6;
	int N = (Mi+1)*Mo;

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
}

void init_global(int T, int Mi, int Mo){
	gl_init_nthreads(T);
	double pe = 1.53;
	gl_init_grid(Mi, Mo, 100, 100, pe);
	gl_init_boxsize(2, 2);
	gl_init_Re(100);
}

const double xaxis1 = 30;
const double xaxis2 = 2e4;
const double yaxis1 = 1e-16;
const double yaxis2 = 1e-4;

/*
 * a = parameter  in (D*D-a*a)u = f
 * flag = SINPIY: u = sin(pi*y)
 *      = FHYY:   u = sin(pi*y)/(1+400y*y)
 *      = ONE:    f = 1
 *      = ... see above for more
 */
void graphsinglegrid(double a, enum bvp_ex flag){
	int Mo = 1;
	int Milist[] = {32, 64, 128, 
			192, 256, 512, 
			1024, 1024*4, 1024*16};
	const int count = 9;
	double mlist[count], 
		e1[count], e2[count], e3[count], e4[count], e5[count];
	
	for(int i=0; i < count; i++)
		mlist[i] = Milist[i];
	
	for(int i=0; i < count; i++){
		int Mi = Milist[i];
		init_global(1, Mi, Mo);
		BVPErrors bvpe(a, Mi, Mo);
		bvpe.inity(gl_y);
		bvpe.initexact(flag);
		bvpe.initcheb(gl_oy, gl_chebi);
		bvpe.computesolns();
		bvpe.computeerrors();
		e1[i] = bvpe.error_u();
		e2[i] = bvpe.error_du1();
		e3[i] = bvpe.error_du2();
		e4[i] = bvpe.error_ddu1();
		e5[i] = bvpe.error_ddu2();
		gl_exit();
	}

	//array_show(mlist, count, "mlist");
	//array_show(e1, count, "elist");
	PyPlot plt("singlegrid");
	plt.loglog(mlist, e1, count);
	plt.linewidth("2");

	plt.loglog(mlist, e2, count);
	plt.linestyle("--");
	plt.linewidth("2");

	plt.loglog(mlist, e3, count);
	plt.linestyle("--");
	plt.linewidth("2");
	plt.markertype("o");
	plt.markercolor("black");
	
	/*
	plt.loglog(mlist, e4, count);
	plt.linestyle("-.");
	plt.linewidth("2");

	plt.loglog(mlist, e5, count);
	plt.linestyle("-.");
	plt.linewidth("2");
	plt.markertype("o");
	plt.markercolor("black");
	*/
	
	plt.axis(xaxis1, xaxis2, yaxis1, yaxis2);
	plt.title(l_title);
	plt.show();
	//plt.savescript();
}

/*
 * a = parameter  in (D*D-a*a)u = f
 * flag = SINPIY: u = sin(pi*y)
 *      = FHYY:   u = sin(pi*y)/(1+400y*y)
 *      = ONE:    f = 1
 */
void graphfixedN(double a, enum bvp_ex flag){
	int Milist[] = {5, 8, 16, 32, 
			64, 128, 192, 256, 
			512, 1024, 1024*2, 1024*4, 
			1024*8, 9999, 1024*16};
	const int count = 15;
	double mlist[count], 
		e1[count], e2[count], e3[count], e4[count], e5[count];
	for(int i=0; i < count; i++)
		mlist[i] = Milist[i];
	for(int i=0; i < count; i++){
		int Mi = Milist[i];
		int Mo = 20000/(Mi+1);
		init_global(1, Mi, Mo);
		BVPErrors bvpe(a, Mi, Mo);
		bvpe.inity(gl_y);
		bvpe.initexact(flag);
		bvpe.initcheb(gl_oy, gl_chebi);
		bvpe.computesolns();
		bvpe.computeerrors();
		e1[i] = bvpe.error_u();
		e2[i] = bvpe.error_du1();
		e3[i] = bvpe.error_du2();
		e4[i] = bvpe.error_ddu1();
		e5[i] = bvpe.error_ddu2();
		gl_exit();
	}

	//array_show(mlist, count, "mlist");
	//array_show(e1, count, "elist");
	PyPlot plt("fixedN10000");
	plt.loglog(mlist, e1, count);
	plt.linewidth("2");

	plt.loglog(mlist, e2, count);
	plt.linestyle("--");
	plt.linewidth("2");

	plt.loglog(mlist, e3, count);
	plt.linestyle("--");
	plt.linewidth("2");
	plt.markertype("o");
	plt.markercolor("black");
	
	/*
	plt.loglog(mlist, e4, count);
	plt.linestyle("-.");
	plt.linewidth("2");

	plt.loglog(mlist, e5, count);
	plt.linestyle("-.");
	plt.linewidth("2");
	plt.markertype("o");
	plt.markercolor("black");
	*/
	
	plt.axis(5, xaxis2, yaxis1, yaxis2);
	plt.title(l_title);
	plt.show();
	//plt.savescript();
}

/*
 * flag = SINPIY: u = sin(pi*y)
 *      = FHYY:   u = sin(pi*y)/(1+400*y*y)
 *      = ONE:    f = 1 
 */
void graphfixedMi(double a, int Mi, enum bvp_ex flag){
	assrt(Mi<300);
	int Nlist[] = {100, 400, 500, 
		       1000, 2000, 5000, 
		       8000, 10000, 20000};
	const int count = 9;
	double nlist[count], 
		e1[count], e2[count], e3[count], e4[count], e5[count];
	for(int i=0; i < count; i++)
		nlist[i] = Nlist[i];
	for(int i=0; i < count; i++){
		int Mo = Nlist[i]/Mi;
		init_global(1, Mi, Mo);
		BVPErrors bvpe(a, Mi, Mo);
		bvpe.inity(gl_y);
		bvpe.initexact(flag);
		bvpe.initcheb(gl_oy, gl_chebi);
		bvpe.computesolns();
		bvpe.computeerrors();
		e1[i] = bvpe.error_u();
		e2[i] = bvpe.error_du1();
		e3[i] = bvpe.error_du2();
		e4[i] = bvpe.error_ddu1();
		e5[i] = bvpe.error_ddu2();
		gl_exit();
	}

	//array_show(nlist, count, "nlist");
	//array_show(e1, count, "elist");
	PyPlot plt("fixedN10000");
	plt.loglog(nlist, e1, count);
	plt.linewidth("2");

	plt.loglog(nlist, e2, count);
	plt.linestyle("--");
	plt.linewidth("2");

	plt.loglog(nlist, e3, count);
	plt.linestyle("--");
	plt.linewidth("2");
	plt.markertype("o");
	plt.markercolor("black");
	
	/*
	plt.loglog(nlist, e4, count);
	plt.linestyle("-.");
	plt.linewidth("2");

	plt.loglog(nlist, e5, count);
	plt.linestyle("-.");
	plt.linewidth("2");
	plt.markertype("o");
	plt.markercolor("black");
	*/

	plt.axis(xaxis1*2, xaxis2*1.1, yaxis1, yaxis2);
	plt.title(l_title);
	plt.show();
}

/*
 * differentiation errors depend greatly on the function being differentiated
 * (D*D-a*a)u = f
 * 1. sinpiy, a = 10, single grid
 * 2. sinkpiy, a = 10, single grid
 */
void show_diff_errors(){
	double a[2] = {10.0, 10.0};
	enum bvp_ex flag[2] = {SINPIY, SINKPIY};

	l_title = "diff errors depend upon function-->sin(pi*y)/single grid";
	graphsinglegrid(a[0], flag[0]);

	l_title = "diff errors depend upon function-->sin(k*pi*y)/single grid";
	graphsinglegrid(a[1], flag[1]);
}

/*
 * differentiation errors in the singularly perturbed setting
 * (D*D-a*a)u = f
 * 1. sinpiy, a = 10^6, single grid
 * 2. sinkpiy, a = 10^6, single grid
 */
void show_diff_largea_errors(){
	double a[2] = {2e4, 2e4};
	enum bvp_ex flag[2] = {SINPIY, SINKPIY};

	l_title = "a = 2e4, sin(pi*y)/single grid";
	graphsinglegrid(a[0], flag[0]);

	l_title = "a = 2e4, sin(k*pi*y)/single grid";
	graphsinglegrid(a[1], flag[1]);
}

/*
 * (D*D-a*a)u = f
 * when spectral integration is based on expanding du (or ddu)
 * the derv (du or ddu) is accurate as u even with large number of cheb points.
 * the two plots here show that this phenomenon exists and that it is partly 
 * a mirage.
 * 1. sinpiy, a=10,   Mi=32 , increasing Mo
 * 2. sinpiy, a=1e6, ditto
 * 1: shows the phenomenon with fixed Mi
 * 2: shows the phenomenon recedes when a is large
 */
void show_accurate_diff_artifact(){
	double a[2]         = {  10.0,  1e6};
	enum bvp_ex flag[2] = {  SINPIY,  SINPIY}; 

	int Mi = 32;
	const char *title[2]  = {"a=10, Mi = 32 fixed, sin(pi*y)",
				   "a=1e6, Mi = 32 fixed, sin(pi*y)"};
	for(int i=0; i < 2; i++){
		l_title = title[i];
		graphfixedMi(a[i], Mi, flag[i]);
	}
}


/*
 * show that a smaller inner grid can be  better
 */
void show_innergriduse(){
	l_title = "smllr innr grd can be bttr, fixed N, a=1e6, sin(pi*y)";
	graphfixedN(1e6, SINPIY);
	
	l_title = "smllr innr grd can be bttr, fixed N, a=1e6, sin(k*pi*y)";
	graphfixedN(1e6, SINKPIY);
}


/*
 * show that if there are only a few grid points per wavelength
 * it may be better to use a larger num of cheb points
 */
void show_gridptperwave(){
	l_title = "few grd pts per wvln,a=1e6,sin(k*piy),fixed N,k=N/6";
	graphfixedN(1e6, SINKPIY);
}

int  main(){
	//show_diff_errors();
	show_diff_largea_errors();
	//show_accurate_diff_artifact();
	//show_innergriduse();
	//show_gridptperwave();
	//l_title = "a=1e5, nsmth, fixedN";
	//graphfixedN(1e5, NSMTH);
}
