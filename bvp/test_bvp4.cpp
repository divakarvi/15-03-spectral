#include <iostream>
#include <cmath>
#include "../utils/utils.hh"
#include "../utils/Table.hh"
#include "../pyplot/pyplot.hh"
#include "bvp4.hh"

void test_bvp4fac(){
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

int main(){
	test_bvp4fac();
}
