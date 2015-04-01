#include <iostream>
#include <cmath>
#include "../utils/utils.hh"
#include "../utils/Table.hh"
#include "../pyplot/pyplot.hh"
#include "bvpi.hh"
#include "bvp4.hh"

void bvp4facA(){
	int M = 32;
	double alpha = 1e3;
	double beta = 1e6;
	
	BVP4fac bvp(alpha, beta, M);

	double f[M+1];
	double u[M+1];
	double ue[M+1];
	double chebi[M+1];

	cheb(chebi, M);
	for(int j=0; j <= M; j++){
		double x = chebi[j];
		ue[j] = sin(PI*x)*sin(PI*x);
		f[j] = -8*PI*PI*PI*PI*cos(2*PI*x) 
			- (alpha*alpha + beta*beta)*(2*PI*PI*cos(2*PI*x))
			+ alpha*alpha*beta*beta*ue[j];
	}
	
	bvp.solve(f, u);
	array_show(u+M-9, 10, "u last ten entries");
	array_show(ue+M-9, 10, "ue last ten entries");

	{
		PyPlot plt("uue");
		plt.plot(chebi, u, M+1);
		plt.linestyle("None");
		plt.markertype("o");
		plt.markercolor("black");
		
		plt.plot(chebi, ue, M+1);
		plt.title("uue");
		plt.show();
	}

	array_diff(u, ue, M+1);
	double emax = array_max(u, M+1);
	
	std::cout<<"BVP4fac test: error in sin(x) = "<<emax<<std::endl;
}

void bvp4facB(){
	int M = 32;
	double alpha = 1.2;
	double beta = 2.0;
	
	BVP4fac bvp(alpha, beta, M);

	double f[M+1];
	double u[M+1];
	double ue[M+1];
	double chebi[M+1];

	cheb(chebi, M);
	for(int j=0; j <= M; j++){
		double a = alpha, b = beta;
		double x = chebi[j];
		ue[j] = -7.0*x*(5.0*x*x - 3.0)/9 
			+ 5.0*x*(63.0*x*x*x*x - 70.0*x*x + 15)/72 
			+ x;
		f[j] = 35.0*x*(
			       pow(a,2)*pow(b,2)*(pow(x,4) - 2.0*pow(x,2) + 1) 
			       - 4.0*(a*a + b*b)*(5*x*x - 3) 
			       + 120.0)/8.0;
	}
	
	bvp.solve(f, u);
	array_show(u+M-9, 10, "u last ten entries");
	array_show(ue+M-9, 10, "ue last ten entries");

	{
		PyPlot plt("uue");
		plt.plot(chebi, u, M+1);
		plt.linestyle("None");
		plt.markertype("o");
		plt.markercolor("black");
		
		plt.plot(chebi, ue, M+1);
		plt.title("uue");
		plt.show();
	}

	array_diff(u, ue, M+1);
	double emax = array_max(u, M+1);
	
	std::cout<<"BVP4fac test: error for odd fn = "<<emax<<std::endl;
}

void bvp4si_solvep(){
	double a = 1.0;
	double b = 2.0;
	int M = 64;

	BVP4si bvp(-a*a-b*b, a*a*b*b, M);

	double f[M+1], ue[M+1], u[M+1], xx[M+1];
	cheb(xx, M);
	for(int j=0; j <= M; j++)
		f[j] = ue[j] = u[j] = 0.0;
	/* 
	for(int j=0; j <= M; j++){//[u = T5] 
		double x = xx[j];
		ue[j] = 16*pow(x,5) - 20*pow(x,3) + 5*x;
		f[j] = a*a*b*b*(16*pow(x,5) 
				- 20*pow(x,3) + 5*x) 
			- 40*x*(a*a + b*b)*(8*x*x - 3) + 1920*x;
	}
	*/
	
	
	for(int k=4; k <= 63; k++)
		for(int j=0; j <= M; j++){
			double t = j*PI/M;
			ue[j] += k*cos(k*t);

			f[j] += a*a*b*b*k*cos(k*t);
			
			if(j != 0 && j != M)
				f[j] -= (a*a+b*b)*
					k*
					k*(-k*sin(t)*cos(k*t) + 
					   sin(k*t)*cos(t))/pow(sin(t),3);
			else{ 
				int sign = (k%2 == 1 && j == M)?-1:1;
				f[j] -= (a*a+b*b)*
					k*
					sign*k*k*(k*k - 1.0)/3;
			}
			
			
			if(j != 0 && j != M)
				f[j] += k*
					k*(k*k*k*cos(k*t) 
					   - 6*k*k*sin(k*t)/tan(t) 
					   + 11*k*cos(k*t) 
					   - 15*k*cos(k*t)/pow(sin(t),2) 
					   + 9*sin(k*t)/tan(t) 
					   + 15*sin(k*t)/pow(tan(t),3))
					/pow(sin(t),4);
			else{
				int sign = (k%2 == 1 && j == M)?-1:1;
				f[j] += 1.0*sign*k*
					k*k*
					(pow(k*1.0,6) 
					 - 14*pow(k*1.0, 4) 
					 + 49*k*k - 36.0)/105.0;
				}
			}

	double dup1, dum1;
	bvp.solvep(f, u, dup1, dum1);
	
	{
		PyPlot plt("bvp4si");
		plt.plot(xx, u, M+1);
		plt.linewidth("3");
		plt.plot(xx, ue, M+1);
		plt.linestyle("None");
		plt.markertype("o");
		plt.markercolor("black");
		plt.show();
	 }

	array_diff(ue, u, M+1);
	printf("error = %e \n", array_max(ue, M+1)/array_max(u, M+1));

	DCT dct(M);
	dct.fwd(u);
	array_show(u, 15, "uf");
}

void bvp4siA(){
	int M = 32;
	double alpha = 1.0e3;
	double beta = 2.0e5;

	BVP4si bvp(-alpha*alpha-beta*beta, alpha*alpha*beta*beta, M);

	double f[M+1];
	double u[M+1];
	double ue[M+1];
	double chebi[M+1];

	cheb(chebi, M);
	for(int j=0; j <= M; j++){
		double x = chebi[j];
		ue[j] = sin(PI*x)*sin(PI*x);
		f[j] = -8*PI*PI*PI*PI*cos(2*PI*x) 
			- (alpha*alpha + beta*beta)*(2*PI*PI*cos(2*PI*x))
			+ alpha*alpha*beta*beta*ue[j];
	}

	bvp.solve(f, u);
	array_show(u+M-9, 10, "u last ten entries");
	array_show(ue+M-9, 10, "ue last ten entries");

	{
		PyPlot plt("uue");
		plt.plot(chebi, u, M+1);
		plt.linestyle("None");
		plt.markertype("o");
		plt.markercolor("black");
		
		plt.plot(chebi, ue, M+1);
		plt.title("uue");
		plt.show();
	 }
	
	array_diff(u, ue, M+1);
	double emax = array_max(u, M+1);
	
	std::cout<<"BVP4si test: error in sin(x) = "<<emax<<std::endl;
}

void bvp4siB(){
	int M = 32;
	double alpha = 1.0e3;
	double beta = 2.0e5;

	BVP4si bvp(-alpha*alpha-beta*beta, alpha*alpha*beta*beta, M);

	double f[M+1];
	double u[M+1];
	double ue[M+1];
	double chebi[M+1];

	cheb(chebi, M);
	for(int j=0; j <= M; j++){
		double a = alpha, b = beta;
		double x = chebi[j];
		ue[j] = -7.0*x*(5.0*x*x - 3.0)/9 
			+ 5.0*x*(63.0*x*x*x*x - 70.0*x*x + 15)/72 
			+ x;
		f[j] = 35.0*x*(
			       pow(a,2)*pow(b,2)*(pow(x,4) - 2.0*pow(x,2) + 1) 
			       - 4.0*(a*a + b*b)*(5*x*x - 3) 
			       + 120.0)/8.0;
	}

	bvp.solve(f, u);
	array_show(u+M-9, 10, "u last ten entries");
	array_show(ue+M-9, 10, "ue last ten entries");

	{
		PyPlot plt("uue");
		plt.plot(chebi, u, M+1);
		plt.linestyle("None");
		plt.markertype("o");
		plt.markercolor("black");
		
		plt.plot(chebi, ue, M+1);
		plt.title("uue");
		plt.show();
	 }
	
	array_diff(u, ue, M+1);
	double emax = array_max(u, M+1);
	
	std::cout<<"BVP4si test: error in odd fn = "<<emax<<std::endl;
}

/*
 * disable Legendre
 */
void bvp4pgA(){
	double a = 2.0;
	double b = 0.0;
	int M =  12;
	BVP4pg bvp(a, b, M);

	double f[M+1];
	for(int i=0; i <= 12; i++)
		f[i] = i;
	double u[M+1];
	bvp.solve(f, u);
	
}

void bvp4pgB(){
	int M = 64;
	double alpha = 1.0;
	double beta = 2.0;

	BVP4pg bvp(alpha*alpha*beta*beta, alpha*alpha+beta*beta, M);

	double f[M+1];
	double u[M+1];
	double ue[M+1];
	double chebi[M+1];

	cheb(chebi, M);
	int k = 0;
	double d = 1.0/sqrt(2.0*(2*k+3)*(2*k+3)*(2*k+5));
	printf("coeffs = %e, %e, %e\n", d, -d*10/7, d*3/7);
	for(int j=0; j <= M; j++){
		double a = alpha;
		double b = beta;
		double x = chebi[j];
		ue[j] = 15.0*pow(x,4)/8 - 15.0*pow(x,2)/4 + 15.0/8;
		ue[j] *= d;
		f[j] = 15.0*pow(a,2)*pow(b,2)*
			(pow(x,4) - 2.0*pow(x,2) + 1.0)/8 
			- 15.0*(pow(a,2) + pow(b,2))*(3.0*pow(x,2) - 1.0)/2 
			+ 45.0;
		f[j] = d*f[j];
	}

	bvp.solve(f, u);

	{
		PyPlot plt("uue");
		plt.plot(chebi, u, M+1);
		plt.linestyle("None");
		plt.markertype("o");
		plt.markercolor("black");
		
		plt.plot(chebi, ue, M+1);
		plt.title("uue");
		plt.show();
	 }

	array_diff(u, ue, M+1);
	double emax = array_max(u, M+1);
	
	std::cout<<"BVP4pg test: error in even soln psi_0 = "<<emax<<std::endl;
}

void bvp4pgC(){
	int M = 64;
	double alpha = 1.0e3;
	double beta = 2.0e5;

	BVP4pg bvp(alpha*alpha*beta*beta, alpha*alpha+beta*beta, M);

	double f[M+1];
	double u[M+1];
	double ue[M+1];
	double chebi[M+1];

	cheb(chebi, M);
	for(int j=0; j <= M; j++){
		double x = chebi[j];
		ue[j] = sin(PI*x)*sin(PI*x);
		f[j] = -8*PI*PI*PI*PI*cos(2*PI*x) 
			- (alpha*alpha + beta*beta)*(2*PI*PI*cos(2*PI*x))
			+ alpha*alpha*beta*beta*ue[j];
	}

	bvp.solve(f, u);

	{
		PyPlot plt("uue");
		plt.plot(chebi, u, M+1);
		plt.linestyle("None");
		plt.markertype("o");
		plt.markercolor("black");
		
		plt.plot(chebi, ue, M+1);
		plt.title("uue");
		plt.show();
	 }

	array_diff(u, ue, M+1);
	double emax = array_max(u, M+1);
	
	std::cout<<"BVP4pg test: error in sin(x) = "<<emax<<std::endl;
}

void bvp4pgD(){
	int M = 64;
	double alpha = 1.0e3;
	double beta = 2.0e5;

	BVP4pg bvp(alpha*alpha*beta*beta, alpha*alpha+beta*beta, M);

	double f[M+1];
	double u[M+1];
	double ue[M+1];
	double chebi[M+1];

	cheb(chebi, M);
	for(int j=0; j <= M; j++){
		double a = alpha, b = beta;
		double x = chebi[j];
		ue[j] = -7.0*x*(5.0*x*x - 3.0)/9 
			+ 5.0*x*(63.0*x*x*x*x - 70.0*x*x + 15)/72 
			+ x;
		f[j] = 35.0*x*(
			       pow(a,2)*pow(b,2)*(pow(x,4) - 2.0*pow(x,2) + 1) 
			       - 4.0*(a*a + b*b)*(5*x*x - 3) 
			       + 120.0)/8.0;
	}

	bvp.solve(f, u);

	{
		PyPlot plt("uue");
		plt.plot(chebi, u, M+1);
		plt.linestyle("None");
		plt.markertype("o");
		plt.markercolor("black");
		
		plt.plot(chebi, ue, M+1);
		plt.title("uue");
		plt.show();
	 }

	array_diff(u, ue, M+1);
	double emax = array_max(u, M+1);
	
	std::cout<<"BVP4pg test: error in odd soln = "<<emax<<std::endl;
}

int main(){
	bvp4facA();
	bvp4facB();
	//bvp4si_solvep();
	//bvp4siA();
	//bvp4siB();
	//bvp4pgA();
	//bvp4pgB();
	//bvp4pgC();
	//bvp4pgD();
}
