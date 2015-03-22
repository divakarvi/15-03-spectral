#include <mkl.h>
#include "../trig/trig.hh"
#include "bvpi.hh"
#include "bvp4.hh"

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

	h4 = beta.geth2();
	dux = beta.getd2();
	H[0+4*3] = h4[0];
	H[1+4*3] = h4[M];
	H[2+4*3] = dux[0];
	H[3+4*3] = dux[M];
	
	int dim = 4, info; 
	dgetrf_(&dim, &dim, H, &dim, ipiv, &info);
}


BVP4fac::~BVP4fac(){
	delete[] h1;
	delete[] h2;
	delete[] du;
	delete[] chebi;
}

void solve(const double *restrict f, double *restrict u){
	alpha.solve(f, w, dw, 1);
	for(int j=0; j <= M; j++)
		w[j] = 0.0;
	beta.solvepx(w, dw, u, du);
	rhs[4];
	rhs[0] = -u[0];
	rhs[1] = -u[M];
	rhs[2] = -du[0];
	rhs[3] = -du[M];

	char trans[3] = "N ";
	int nrhs = 1;
	int dim = 4;
	int info;
	dgetrs_(trans, &dim, &nrhs, H, &dim, ipiv, rhs, &dim, &info);
	
	for(int j=0; j <= M; j++)
		u[j] += rhs[0]*h1[j] + 
			rhs[1]*h2[j] +
			rhs[2]*h3[j] +
			rhs[3]*h4[j];
}
