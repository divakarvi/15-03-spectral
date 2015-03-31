#include <iostream>
#include <cmath>
#include <omp.h>
#include "../utils/utils.hh"
#include "../utils/Table.hh"
#include "../utils/TimeStamp.hh"
#include "../pyplot/pyplot.hh"
#include "bvpi.hh"
#include "bvp4.hh"

#define NTHREADS 16

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


enum TRIAL {INCACHE, OUTCACHE, THREAD};
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

void incache(int count, enum SOLVER solver, int M, struct ab_struct ab, 
	     double *f, double *u){
	switch(solver){
	case FAC:
		{
			BVP4fac bvp4(ab.a, ab.b, M);
			for(long i=0; i < count; i++)
				bvp4.solve(f, u);
		}
		break;
	case SI:
		{
			BVP4si bvp4(ab.a, ab.b, M);
			for(long i=0; i < count; i++)
				bvp4.solve(f, u);
		}
		break;
	case PG:
		{
			BVP4pg bvp4(ab.a, ab.b, M);
			for(long i=0; i < count; i++)
				bvp4.solve(f, u);
		}
		break;
	}
}


void outcache(int count, enum SOLVER solver, int M, struct ab_struct ab, 
	      double *f, double *u){
	switch(solver){
	case FAC:
		{
			BVP4fac bvp4(ab.a, ab.b, M);
			for(long i=0; i < count; i++)
				bvp4.solve(f+i*(M+1), u+i*(M+1));
		}
		break;
	case SI:
		{
			BVP4si bvp4(ab.a, ab.b, M);
			for(long i=0; i < count; i++)
				bvp4.solve(f+i*(M+1), u+i*(M+1));
		}
		break;
	case PG:
		{
			BVP4pg bvp4(ab.a, ab.b, M);
			for(long i=0; i < count; i++)
				bvp4.solve(f+i*(M+1), u+i*(M+1));
		}
		break;
	}
}

void thread(int count, enum SOLVER solver, int M, struct ab_struct ab, 
	      double *f, double *u){
	switch(solver){
	case FAC:
#pragma omp parallel num_threads(NTHREADS)
		{
			int t = omp_get_thread_num();
			BVP4fac bvp4(ab.a, ab.b, M);
			for(long i=count/NTHREADS*t; 
			    i < count/NTHREADS*(t+1); i++)
				bvp4.solve(f+i*(M+1), 
					   u+i*(M+1));
		}
		break;
	case SI:
#pragma omp parallel  num_threads(NTHREADS)
		{
			int t = omp_get_thread_num();
			BVP4si bvp4(ab.a, ab.b, M);
			for(long i=count/NTHREADS*t; 
			    i < count/NTHREADS*(t+1); i++)
				bvp4.solve(f+i*(M+1), 
					   u+i*(M+1));
		}
		break;
	case PG:
#pragma omp parallel  num_threads(NTHREADS)
		{
			int t = omp_get_thread_num();
			BVP4pg bvp4(ab.a, ab.b, M);
			for(long i=count/NTHREADS*t; 
			    i < count/NTHREADS*(t+1); i++)
				bvp4.solve(f+i*(M+1), 
					   u+i*(M+1));
		}
		break;
	}
}

/*
 * return num of cycles per solve/M
 */
double get_cycles(enum TRIAL trial, enum SOLVER solver,
		  double alpha, double beta, int M){
	struct ab_struct ab;
	alphabeta2ab(solver, alpha, beta, &ab);

	long memone = 2*(M+1)*8;
	long count = 1l*1000*1000*1000/memone;
	printf("count = %ld\n", count);
	double *u = (double *)malloc(count*(M+1)*sizeof(double));
	double *f = (double *)malloc(count*(M+1)*sizeof(double));
	double *usave = u, *fsave = f;
	double xx[M+1];
	cheb(xx, M);
	
	/*
	 * initialize f[] and u[]
	 */
	if(trial != THREAD){
		for(long i=0; i < count; i++){
			for(int j=0; j <= M; j++){
				double x = xx[j];
				f[j+i*(M+1)] = f1(x, alpha, beta);
				u[j+i*(M+1)] = u1(x, alpha, beta);
			}
		}
	}
	else{
#pragma omp parallel for			\
	num_threads(NTHREADS)
		for(long i=0; i < count; i++){
			for(int j=0; j <= M; j++){
				double x = xx[j];
				f[j+i*(M+1)] = f1(x, alpha, beta);
				u[j+i*(M+1)] = u1(x, alpha, beta);
			}
		}
	}

	TimeStamp clk;
	clk.tic();
	switch(trial){
	case INCACHE:
		incache(count, solver, M, ab, f, u);
		break;
	case OUTCACHE:
		outcache(count, solver, M, ab, f, u);
		break;
	case THREAD:
		thread(count, solver, M, ab, f, u);
		break;
	}
	double cycles = clk.toc();
	free(usave);
	free(fsave);


	cycles = cycles/count/(M+1);
	printf("M = %d, cycles = %f \n", M, cycles);
	return cycles;
}

void get_cyclelist(enum TRIAL trial, enum SOLVER solver,
		     double alpha, double beta, 
		     int *Mlist, double *cyclelist, int len){
	for(int i=0; i < len; i++)
		cyclelist[i] = get_cycles(trial, solver, alpha, beta,
					  Mlist[i]);
}

int Mlist[11] = {64, 128, 256, 512, 1024, 2048, 4096, 
		 8192, 16384, 32768, 32768*2};

void plot(enum TRIAL trial, double alpha, double beta){
	PyPlot plt("usinpixsq");//dbg here
	double cyclelist[11];
	double mlist[11];
	for(int i=0; i < 11; i++)
		mlist[i] = Mlist[i];

	get_cyclelist(trial, FAC, alpha, beta, Mlist, cyclelist, 11);
	plt.semilogx(mlist, cyclelist, 11);
	plt.linestyle("-");
	plt.linecolor("black");
	plt.linewidth("3");

	get_cyclelist(trial, SI, alpha, beta, Mlist, cyclelist, 11);
	plt.semilogx(mlist, cyclelist, 11);
	plt.linestyle("--");
	plt.linecolor("black");
	plt.linewidth("3");
	
	if (trial != THREAD) {
		get_cyclelist(trial, PG, alpha, beta, Mlist, cyclelist, 11);
		plt.semilogx(mlist, cyclelist, 11);
		plt.linestyle(":");
		plt.linecolor("black");
		plt.linewidth("3");
	}

	plt.show();
}

int main(){
	kmp_set_defaults("KMP_AFFINITY=compact");
	//plot(INCACHE, 1.0, 2.0);
	plot(OUTCACHE, 1.0, 2.0);
	plot(THREAD, 1.0, 2.0);
}
