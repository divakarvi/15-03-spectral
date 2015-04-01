#include <iostream>
#include <cmath>
#include "../utils/utils.hh"
#include "../utils/Table.hh"
#include "../pyplot/pyplot.hh"
#include "bvpi.hh"
#include "bvp4.hh"


/*
 * (D^2 - alpha^2)(D^2-beta^2) u = f
 * f1, u1: u = sin(pi*x)
 * f2, u2: f = 1*alpha^2*beta^2
 */
double f1(double x, double alpha, double beta){
	double ans  = 
		-8*PI*PI*PI*PI*cos(2*PI*x) 
		- (alpha*alpha + beta*beta)*(2*PI*PI*cos(2*PI*x))
		+ alpha*alpha*beta*beta*sin(PI*x)*sin(PI*x);
	return ans;
}

double u1(double x, double alpha, double beta){
	double ans = sin(PI*x)*sin(PI*x);
	return ans;
}

double f2(double x, double alpha, double beta){
	double ans  = alpha*alpha*beta*beta;
	return ans;
}

double u2(double x, double alpha, double beta){
	double ans;
	double y = x;
	ans = pow(alpha, -0.2e1) * pow(beta, -0.2e1) - (-exp(beta * (-0.1e1 + y)) + exp(-0.2e1 * alpha - beta + beta * y)) / alpha * pow(beta, -0.2e1) / (-alpha + beta + exp(-0.2e1 * alpha) * alpha + exp(-0.2e1 * alpha) * beta - exp(-0.2e1 * beta) * alpha - exp(-0.2e1 * beta) * beta + exp(-0.2e1 * alpha - 0.2e1 * beta) * alpha - exp(-0.2e1 * alpha - 0.2e1 * beta) * beta) - (exp(alpha * (-0.1e1 + y)) - exp(-alpha - 0.2e1 * beta + alpha * y)) / beta * pow(alpha, -0.2e1) / (-alpha + beta + exp(-0.2e1 * alpha) * alpha + exp(-0.2e1 * alpha) * beta - exp(-0.2e1 * beta) * alpha - exp(-0.2e1 * beta) * beta + exp(-0.2e1 * alpha - 0.2e1 * beta) * alpha - exp(-0.2e1 * alpha - 0.2e1 * beta) * beta) - (-exp(-beta * (0.1e1 + y)) + exp(-0.2e1 * alpha - beta - beta * y)) / alpha * pow(beta, -0.2e1) / (-alpha + beta + exp(-0.2e1 * alpha) * alpha + exp(-0.2e1 * alpha) * beta - exp(-0.2e1 * beta) * alpha - exp(-0.2e1 * beta) * beta + exp(-0.2e1 * alpha - 0.2e1 * beta) * alpha - exp(-0.2e1 * alpha - 0.2e1 * beta) * beta) - (exp(-alpha * (0.1e1 + y)) - exp(-alpha - 0.2e1 * beta - alpha * y)) / beta * pow(alpha, -0.2e1) / (-alpha + beta + exp(-0.2e1 * alpha) * alpha + exp(-0.2e1 * alpha) * beta - exp(-0.2e1 * beta) * alpha - exp(-0.2e1 * beta) * beta + exp(-0.2e1 * alpha - 0.2e1 * beta) * alpha - exp(-0.2e1 * alpha - 0.2e1 * beta) * beta);
	ans = alpha*alpha*beta*beta*ans;
	return ans;
}


enum TRIAL {SINPIX, ONE};
enum SOLVER {FAC, SI, PG}; /*factorized bvp4, spec intgn, petrov-galerkin*/

struct ab_struct{
	double a;
	double b;
};

void alphabeta2ab(enum  SOLVER solver, double alpha, double beta,
		 struct ab_struct *ab){
	switch(solver){
	case FAC:
		ab->a = alpha;
		ab->b = beta;
		break;
	case SI:
		ab->a = -alpha*alpha - beta*beta;
		ab->b = alpha*alpha*beta*beta;
		break;
	case PG:
		ab->a = alpha*alpha*beta*beta;
		ab->b = alpha*alpha + beta*beta;
		break;
	}
}

double find_error(enum TRIAL trial, enum SOLVER solver,
		  double alpha, double beta, int M){
	struct ab_struct ab;
	alphabeta2ab(solver, alpha, beta, &ab);

	double xx[M+1];
	cheb(xx, M);
	double u[M+1], f[M+1], ue[M+1];
	for(int j=0; j <= M; j++){
		double x = xx[j];
		switch(trial){
		case SINPIX:
			f[j] = f1(x, alpha, beta);
			ue[j] = u1(x, alpha, beta);
			break;
		case ONE:
			f[j] = f2(x, alpha, beta);
			ue[j] = u2(x, alpha, beta);
			break;
		}
		u[j] = 0.0;
	}

	switch(solver){
	case FAC:
		{
			BVP4fac bvp4fac(ab.a, ab.b, M);
			bvp4fac.solve(f, u);
		}
		break;
	case SI:
		{
			BVP4si bvp4si(ab.a, ab.b, M);
			bvp4si.solve(f, u);
		}
		break;
	case PG:
		{
			BVP4pg bvp4pg(ab.a, ab.b, M);
			bvp4pg.solve(f, u);
		}
		break;
	}
	array_diff(ue, u, M+1);
	double err = array_max(ue, M+1);
	return err;
}

void find_errlist(enum TRIAL trial, enum SOLVER solver,
		  double alpha, double beta, 
		  int *Mlist, double *errlist, int len){
	for(int i=0; i < len; i++)
		errlist[i] = find_error(trial, solver, alpha, beta, Mlist[i]);
}

int Mlist[11] = {64, 128, 256, 512, 1024, 2048, 4096, 
		 8192, 16384, 32768, 65536};

void plot(enum TRIAL trial, double alpha, double beta){
	char name[100];
	if (trial == SINPIX)
		sprintf(name, "accuracy_sinpix");
	else if(trial == ONE)
		sprintf(name, "accuracy_blayer");
	//PyPlot plt(name);
	PyPlot plt(name, PLTOFF);//PLTOFF for pdf
	double errlist[11];
	double mlist[11];
	for(int i=0; i < 11; i++)
		mlist[i] = Mlist[i];
	
	find_errlist(trial, FAC, alpha, beta, Mlist, errlist, 11);
	plt.loglog(mlist, errlist, 11);
	plt.linestyle("-");
	plt.linecolor("black");
	plt.linewidth("3");
	
	find_errlist(trial, SI, alpha, beta, Mlist, errlist, 11);
	plt.loglog(mlist, errlist, 11);
	plt.linestyle("--");
	plt.linecolor("black");
	plt.linewidth("3");

	find_errlist(trial, PG, alpha, beta, Mlist, errlist, 11);
	plt.loglog(mlist, errlist, 11);
	plt.linestyle(":");
	plt.linecolor("black");
	plt.linewidth("3");
	
	const char* cmds = "plt.legend(['FAC', 'SI', 'PG'], loc = 'lower left')"
		"\n"
		"ax.set_xticks([2**6, 2**7, 2**8, 2**9, 2**10,"
		"2**11, 2**12, 2**13, 2**14, 2**15, 2**16])"
		"\n"
		"ax.set_xticklabels([r'$2^6$', r'$2^7$', "
		"r'$2^8$', r'$2^9$', r'$2^{10}$', r'$2^{11}$', r'$2^{12}$'," 
		"r'$2^{13}$', r'$2^{14}$', r'$2^{15}$', r'$2^{16}$'])"
		"\n"
		"plt.xlabel('M', fontsize=20)"
		"\n"
		"plt.ylabel('Error', fontsize=14)"
		"\n";
	plt.pycmd(cmds);

	plt.show();
}

int main(){
	plot(SINPIX, 1e3, 1e6);
	plot(ONE, 1e3, 1e6);
}
