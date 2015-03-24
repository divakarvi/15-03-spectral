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
}


BVP4si::~BVP4si(){
	delete[] space;
}

void BVP4si::solvep(const double *restrict ff, double *restrict u, 
		    double& dup1, double& dum1){
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
		dum1 = u[j]*j*j - u[j+1]*(j+1)*(j+1);

	dct.bwd(u);
}
