#include <cmath>
#include <mkl.h>
#include "../utils/utils.hh"
#include "../banded/banded.hh"
#include "bvp4.hh"

#define dk(k) (1.0/sqrt((2.0*(2*(k)+3.0)*(2*(k)+3.0)*(2*(k)+5.0))))
#define ek(k) (2.0/(2*(k)+1.0))
#define gk(k) ((2*(k)+3.0)/(2*(k)+7.0))
#define hk(k) (-(1.0+gk(k)))
#define bk(k)  dk(k)*dk(k)*(ek(k) + hk(k)*hk(k)*ek((k)+2) + gk(k)*gk(k)*ek((k)+4))
#define bkp2(k) dk(k)*dk((k)+2)*(hk(k)*ek((k)+2) + gk(k)*hk((k)+2)*ek((k)+4))
#define bkp4(k) dk(k)*dk((k)+4)*gk((k))*ek((k)+4)
#define ck(k) (-2.0*(2*(k)+3.0)*dk(k)*dk(k)*hk(k))
#define ckp2(k) (-2.0*(2.0*(k)+3)*dk(k)*dk((k)+2))
#define ggk(k) (-2.0*(2*(k)+5)/(2*(k)+7))
#define lk(k) (2.0/(2.0*(k)+1))

BVP4pg::BVP4pg(double aa, double bb, int MM)
	:a(aa), b(bb), M(MM),
	Me((M-3+(1-M%2))/2), Mo((M-3-(1-M%2))/2),
	even(Me), odd(Mo),
	 flt(M) 
{
	space = new double[2*(M+1)];
	
	double *restrict u2 = space;
	for(int k=0; 2*k+4 <= M-4; k++)
		u2[k] = a*bkp4(2*k);
	even.setu2(u2);

	
	double *restrict u1 = space;
	for(int k=0; 2*k+2 <= M-4; k++)
		u1[k] = a*bkp2(2*k) + b*ckp2(2*k);
	even.setu1(u1);

	double *restrict d = space;
	for(int k=0; 2*k <= M-4; k++)
		d[k] = a*bk(2*k) + b*ck(2*k) + 1.0;
	even.setdiag(d);


	even.factor();
	
	u2 = space;
	for(int k=0; 2*k+5 <= M-4; k++)
		u2[k] = a*bkp4(2*k+1);
	odd.setu2(u2);
	
	u1 = space;
	for(int k=0; 2*k+3 <= M-4; k++)
		u1[k] = a*bkp2(2*k+1) + b*ckp2(2*k+1);
	odd.setu1(u1);
	
	d = space;
	for(int k=0; 2*k+1 <= M-4; k++)
		d[k] = a*bk(2*k+1) + b*ck(2*k+1) + 1.0;
	odd.setdiag(d);
	
	odd.factor();
}

BVP4pg::~BVP4pg(){
	delete[] space;
}
 

void BVP4pg::solve(const double *restrict fi, double *restrict u){
	double *restrict f = space;
	for(int j=0; j <= M; j++)
		f[j] = fi[j];
	
	/*
	 * take Legendre transform of f
	 */
	flt.fwd(f);

	double *restrict re = f + (M + 1);
	double *restrict ro = re + (M-3+(1-M%2))/2;
	
	for(int k=0; 2*k <= M-4; k++)
		re[k] = dk(2*k)*(f[2*k]*lk(2*k)
				 + ggk(2*k)*f[2*k+2]*lk(2*k+2)
				 + gk(2*k)*f[2*k+4]*lk(2*k+4));
	even.solve(re, 1);

	for(int k=0; 2*k+1 <= M-4; k++)
		ro[k] = dk(2*k+1)*(f[2*k+1]*lk(2*k+1)
				   + ggk(2*k+1)*f[2*k+1+2]*lk(2*k+3)
				   + gk(2*k+1)*f[2*k+1+4]*lk(2*k+5));
	odd.solve(ro, 1);
	
	int j, k;
	for(k=0; 2*k <= M-4; k++){
		j = 2*k;
		f[j] = dk(j)*re[k]; 
	}
	
	for(k=1; 2*k <= M-4; k++){
		j = 2*k;
		f[j] += dk(j-2)*ggk(j-2)*re[k-1]; 
	}
	for(k=2; 2*k <= M-4; k++){
		j = 2*k;
		f[j] += dk(j-4)*gk(j-4)*re[k-2];
	}

	j = 2*k;
	if (j <= M){
		f[j] = dk(j-2)*ggk(j-2)*re[k-1] 
			+  dk(j-4)*gk(j-4)*re[k-2];
	}
	k += 1;
	j = 2*k;
	if (j <= M){
		f[j] = dk(j-4)*gk(j-4)*re[k-2];
	}

	for(k=0; 2*k+1 <= M-4; k++){
		j = 2*k+1;
		f[j] = dk(j)*ro[k]; 
	}
	for(k=1; 2*k+1 <= M-4; k++){
		j = 2*k + 1;
		f[j] += dk(j-2)*ggk(j-2)*ro[k-1]; 
	}
	for(k=2; 2*k+1 <= M-4; k++){
		j = 2*k + 1;
		f[j] += dk(j-4)*gk(j-4)*ro[k-2];
	}

	j = 2*k+1;
	if (j <= M){
		f[j] = dk(j-2)*ggk(j-2)*ro[k-1] 
			+  dk(j-4)*gk(j-4)*ro[k-2];
	}
	k += 1;
	j = 2*k+1;
	if (j <= M){
		f[j] = dk(j-4)*gk(j-4)*ro[k-2];
	}

	/*
	 * inverse Legendre transform f to u
	 */
	flt.bwd(f);
	for(int i=0; i <= M; i++)
		u[i] = f[i];
}



