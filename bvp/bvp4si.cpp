#include <cmath>
#include <mkl.h>
#include "../utils/utils.hh"
#include "../banded/banded.hh"
#include "../fft/trig.hh"
#include "bvp4.hh"

#define I2m2(n) (1.0/(4.0*(n)*((n)-1.0))) 
#define I2(n) (-1.0/(2.0*((n)-1.0)*((n)+1.0)))
#define I2p2(n) (1.0/(4.0*(n)*((n) + 1.0)))

#define I4m4(n) (1.0/(16.0*(n)*((n) - 3.0)*((n) - 2.0)*((n) - 1.0)))
#define I4m2(n) (-1.0/(4.0*(n)*((n) - 3.0)*((n) - 1.0)*((n) + 1.0)))
#define I4(n) (3.0/(8.0*((n) - 2.0)*((n) - 1.0)*((n) + 1.0)*((n) + 2.0)))
#define I4p2(n) (-1.0/(4.0*(n)*((n) - 1.0)*((n) + 1.0)*((n) + 3.0)))
#define I4p4(n) (1.0/(16.0*(n)*((n) + 1.0)*((n) + 2.0)*((n) + 3.0)))

BVP4si::BVP4si(double aa, double bb, int MM)
	:a(aa), b(bb), M(MM),
	Me((M-4+M%2)/2), Mo((M-4-M%2)/2),
	even(Me), odd(Mo),
	dct(M)
{
	space = new double[2*(M+1)];
	
	double *restrict u2 = space;
	for(int k=0; 2*k+8 < M; k++)
		u2[k] = b*I4p4(2*k+4);
	even.setu2(u2);

	double *restrict u1 = space;
	for(int k=0; 2*k+6 < M; k++)
		u1[k] = a*I2p2(2*k+4) + b*I4p2(2*k+4);
	even.setu1(u1);

	double *restrict d = space;
	for(int k=0; 2*k+4 < M; k++)
		d[k] = 1 + a*I2(2*k+4) + b*I4(2*k+4);
	even.setdiag(d);

	double *restrict l1 = space;
	for(int k=0; 2*k+6 < M; k++)
		l1[k] = a*I2m2(2*k+6) + b*I4m2(2*k+6);
	even.setl1(l1);

	double *restrict l2 = space;
	for(int k=0; 2*k+8 < M; k++)
		l2[k] = b*I4m4(2*k+8);
	even.setl2(l2);

	even.factor();

	u2 = space;
	for(int k=0; 2*k+9 < M; k++)
		u2[k] = b*I4p4(2*k+5);
	odd.setu2(u2);

	u1 = space;
	for(int k=0; 2*k+7 < M; k++)
		u1[k] = a*I2p2(2*k+5) + b*I4p2(2*k+5);
	odd.setu1(u1);

	d = space;
	for(int k=0; 2*k+5 < M; k++)
		d[k] = 1 + a*I2(2*k+5) + b*I4(2*k+5);
	odd.setdiag(d);

	l1 = space;
	for(int k=0; 2*k+7 < M; k++)
		l1[k] = a*I2m2(2*k+7) + b*I4m2(2*k+7);
	odd.setl1(l1);

	l2 = space;
	for(int k=0; 2*k+9 < M; k++)
		l2[k] = b*I4m4(2*k+9);
	odd.setl2(l2);

	odd.factor();

	h1 = new double[M+1];
	h2 = new double[M+1];
	h3 = new double[M+1];
	h4 = new double[M+1];
	double f[M+1], dup1, dum1;
	
	for(int j=0; j <= M; j++)
		f[j] = -b/2;
	this->solvep(f, h1, dup1, dum1, DCTOFF);
	dup1 += 0.0;
	dum1 += 0.0;
	h1[0] = 1.0;
	dct.bwd(h1);
	H[0+0*4] = h1[0];
	H[1+0*4] = h1[M];
	H[2+0*4] = dup1;
	H[3+0*4] = dum1;

	double x[M+1];
	cheb(x, M);

	for(int j=0; j <= M; j++)
		f[j] = -b*x[j];
	this->solvep(f, h2, dup1, dum1, DCTOFF);
	dup1 += 1.0;
	dum1 += 1.0;
	h2[1] = 1.0;
	dct.bwd(h2);
	H[0+1*4] = h2[0];
	H[1+1*4] = h2[M];
	H[2+1*4] = dup1;
	H[3+1*4] = dum1;

	for(int j=0; j <= M; j++)
		f[j] = -(a*4.0 + b*(2*x[j]*x[j]-1));
	this->solvep(f, h3, dup1, dum1, DCTOFF);
	dup1 += 4*1.0;
	dum1 += -4*1.0;
	h3[2] = 1.0;
	dct.bwd(h3);
	H[0+2*4] = h3[0];
	H[1+2*4] = h3[M];
	H[2+2*4] = dup1;
	H[3+2*4] = dum1;

	for(int j=0; j <= M; j++)
		f[j] = -(a*24*x[j] + b*(4*pow(x[j],3)-3*x[j]));
	this->solvep(f, h4, dup1, dum1, DCTOFF);
	dup1 += 9.0;
	dum1 += 9.0;
	h4[3] = 1.0;
	dct.bwd(h4);
	H[0+3*4] = h4[0];
	H[1+3*4] = h4[M];
	H[2+3*4] = dup1;
	H[3+3*4] = dum1;

	int dim = 4, info; 
	dgetrf_(&dim, &dim, H, &dim, ipiv, &info);
}


BVP4si::~BVP4si(){
	delete[] space;
	delete[] h1;
	delete[] h2;
	delete[] h3;
	delete[] h4;
}

void BVP4si::solvep(const double *restrict ff, double *restrict u, 
		    double& dup1, double& dum1, enum DCTONOFF flag){
	double *restrict f = space;
	double *restrict re = f + (M + 1);
	double *restrict ro = re + (M-4+M%2)/2;

	for(int i=0; i <= M; i++)
		f[i] = ff[i];

	dct.fwd(f);

	int k;
	for(k=0; 2*k+8 < M; k++)
		re[k] = f[2*k] * I4m4(2*k+4) +
			f[2*k+2] * I4m2(2*k+4) +
			f[2*k+4] * I4(2*k+4) +
			f[2*k+6] * I4p2(2*k+4) +
			f[2*k+8] * I4p4(2*k+4);
	if(2*k+6 < M)
		re[k] = f[2*k] * I4m4(2*k+4) +
			f[2*k+2] * I4m2(2*k+4) +
			f[2*k+4] * I4(2*k+4) +
			f[2*k+6] * I4p2(2*k+4);
	k = k + 1;
	if(2*k+4 < M)
		re[k] = f[2*k] * I4m4(2*k+4) +
			f[2*k+2] * I4m2(2*k+4) +
			f[2*k+4] * I4(2*k+4);

	even.solve(re, 1);

	for(k=0; 2*k+9 < M; k++)
		ro[k] = f[2*k+1] * I4m4(2*k+5) +
			f[2*k+3] * I4m2(2*k+5) +
			f[2*k+5] * I4(2*k+5) +
			f[2*k+7] * I4p2(2*k+5) +
			f[2*k+9] * I4p4(2*k+5);
	if(2*k+7 < M)
		ro[k] = f[2*k+1] * I4m4(2*k+5) +
			f[2*k+3] * I4m2(2*k+5) +
			f[2*k+5] * I4(2*k+5) +
			f[2*k+7] * I4p2(2*k+5);
	k = k + 1;
	if(2*k+5 < M)
		ro[k] = f[2*k+1] * I4m4(2*k+5) +
			f[2*k+3] * I4m2(2*k+5) +
			f[2*k+5] * I4(2*k+5);

	odd.solve(ro, 1);
		
	u[0] = u[1] = u[2] = u[3] = 0.0;
	for(k=0; 2*k+4 < M; k++)
		u[2*k+4] = re[k];
	for(k=0; 2*k+5 < M; k++)
		u[2*k+5] = ro[k];
	u[M] = 0.0;

	dup1 = 0.0;
	for(int j=1; j < M; j++)
		dup1 += u[j]*j*j;

	dum1 = 0.0;
	for(int j=1; j+1 < M; j=j+2)//OK because u[M] = 0.0
		dum1 += u[j]*j*j - u[j+1]*(j+1)*(j+1);

	if (flag == DCTOFF)
		return;
	dct.bwd(u);
}

void BVP4si::solve(const double *restrict f, double *restrict u){
	double dup1, dum1;
	this->solvep(f, u, dup1, dum1);

	double rhs[4];
	rhs[0] = -u[0];
	rhs[1] = -u[M];
	rhs[2] = -dup1;
	rhs[3] = -dum1;

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
