#include <mkl.h>
#include "../pyplot/pyplot.hh"
#include "../fft/trig.hh"
#include "bvpi.hh"
#include "bvp4.hh"


#undef PLTUPH3

BVP4fac::BVP4fac(double a, double b, int MM)
	:M(MM), alpha(a, M), beta(b, M), dct(M)
{
	chebi = new double[M+1];
	cheb(chebi, M);

	alpha.initcheb(chebi);
	alpha.initdct(dct);
	alpha.inithmg();

	beta.initcheb(chebi);
	beta.initdct(dct);
	beta.inithmg();

	h1 = new double[M+1];
	du = new double[M+1];
	double *f = alpha.geth1();
	beta.solvep(f, h1, du, 1);
	H[0+4*0] = h1[0];
	H[1+4*0] = h1[M];
	H[2+4*0] = du[0];
	H[3+4*0] = du[M];

	h2 = new double[M+1];
	f = alpha.geth2();
	beta.solvep(f, h2, du, 1);
	H[0+4*1] = h2[0];
	H[1+4*1] = h2[M];
	H[2+4*1] = du[0];
	H[3+4*1] = du[M];
		
	h3 = beta.geth1();
	double *dux = beta.getd1();
	H[0+4*2] = h3[0];
	H[1+4*2] = h3[M];
	H[2+4*2] = dux[0];
	H[3+4*2] = dux[M];
	
	{
#ifdef PLTUPH3
		PyPlot plt("h3");
		plt.plot(chebi, h3, M+1);
		plt.title("h3");
		plt.show();
#endif
	}

	h4 = beta.geth2();
	dux = beta.getd2();
	H[0+4*3] = h4[0];
	H[1+4*3] = h4[M];
	H[2+4*3] = dux[0];
	H[3+4*3] = dux[M];
	
	int dim = 4, info; 
	dgetrf_(&dim, &dim, H, &dim, ipiv, &info);

	w = new double[M+1];
}


BVP4fac::~BVP4fac(){
	delete[] h1;
	delete[] h2;
	delete[] du;
	delete[] chebi;
	delete[] w;
}

void BVP4fac::solve(const double *restrict f, double *restrict u){
	alpha.solvep(f, w, du, 1);
	beta.solvep(w, u, du, 1);

	double rhs[4];
	rhs[0] = -u[0];
	rhs[1] = -u[M];
	rhs[2] = -du[0];
	rhs[3] = -du[M];

	char trans[3] = "N ";
	int nrhs = 1;
	int dim = 4;
	int info;
	dgetrs_(trans, &dim, &nrhs, H, &dim, ipiv, rhs, &dim, &info);

	{
#ifdef PLTUPH3
		PyPlot plt("up");
		plt.plot(chebi, u, M+1);
		plt.title("up");
		plt.show();
#endif
	}

	for(int j=0; j <= M; j++)
		u[j] += rhs[0]*h1[j] + 
			rhs[1]*h2[j] +
			rhs[2]*h3[j] +
			rhs[3]*h4[j];
}
