#include "bvpi.hh"
#include "../utils/utils.hh"
#include <cmath>
#include <mkl.h>

void cheb(double *restrict x, int M){
	for(int j=0; j <= M; j++)
		x[j] = cos(j*PI/M);
}

BVPSolvei::BVPSolvei(double aa, int MM)
	:state(0), a(aa),
	 M(MM), Me((M-2+M%2)/2), Mo((M-M%2)/2),
	 even_tri(Me), odd_tri(Mo) 
{
	assrt(state==0);
	state = 1;
	
	assrt(M >= 5);//to ensure Mo, Me >= 2
	
	double tmp[Mo];

	for(int n=2; n <= M-4+M%2; n+=2)
		tmp[n/2-1] = -a*a/(4.0*n*(n+1));
	even_tri.setu1(tmp);

	for(int n=2; n <= M-2+M%2; n+=2)
		tmp[n/2-1] = 1.0 + a*a/(2.0*(n*n-1.0));
	even_tri.setdiag(tmp);

	for(int n=4; n <= M-2+M%2; n+=2)
		tmp[n/2-2] = -a*a/(4.0*n*(n-1));
	even_tri.setl1(tmp);

	even_tri.factor();

	for(int n=1; n <= M-3-M%2; n+=2)
		tmp[(n-1)/2] = -a*a/(4.0*n*(n+1));
	odd_tri.setu1(tmp);

	tmp[0] = 1.0 + a*a/8;
	for(int n=3; n <= M-1-M%2; n+=2)
		tmp[(n-1)/2] = 1.0 + a*a/(2.0*(n*n-1.0));
	odd_tri.setdiag(tmp);

	for(int n=3; n <= M-1-M%2; n+=2)
		tmp[(n-3)/2] = -a*a/(4.0*n*(n-1));
	odd_tri.setl1(tmp);

	odd_tri.factor();

	space = (double *)MKL_malloc(4*(M+1)*sizeof(double), 64);
}

BVPSolvei::~BVPSolvei(){
	MKL_free(space);
}

void BVPSolvei::initcheb(double *chebi){
	assrt(state==1);
	state = 2;
	cheb = chebi;
}

void BVPSolvei::initdct(DCT& dcti){
	assrt(state==2);
	state = 3;
	trig = &dcti;
}

void BVPSolvei::inithmg(){
	assrt(state==3);
	state = 4;
	double f[2*(M+1)];
	
	for(int i=0; i <= M; i++){
		f[i] = a*a/2.0;
		f[M+1+i] = a*a*cheb[i]/2.0;
	}
	
	h1 = space;
	h2 = h1 + (M+1);
	d1 = h2 + (M+1);
	d2 = d1 + (M+1);

	this->solvep(f, h1, d1, 2);

	for(int i=0; i <= M; i++){
		h1[i] += 0.5;
		d2[i] += 0.5;
		h2[i] += 0.5*cheb[i];
	}
}

void BVPSolvei::solvep(const double *restrict f,
		       double *restrict u, double *restrict du, int nrhs){
	assrt(state==4);
	
	/*
	 * copy f to du
	 */
	for(int i=0; i < nrhs*(M+1); i++)
		du[i] = f[i];

	/*
	 * find cheb coeffs of f copied to du
	 */
	for(int i=0; i < nrhs; i++){
		trig->fwd(du+i*(M+1));
		du[M+i*(M+1)] = 0;
	}
	
	{
		double *restrict scre = u;
		double *restrict scro = u + nrhs*Me;
		/*
		 * separate even odd modes of f in du
		 * integrate and pack into scre and scro
		 */
		for(int i=0; i < nrhs; i++){//Mo>=Me
			for(int j = 2; j <= 2*Me; j += 2){
				scre[j/2-1+i*Me]=(du[j-1+i*(M+1)]
						    -du[j+1+i*(M+1)])/(2.0*j);
				scro[j/2-1+i*Mo]=(du[j-2+i*(M+1)]
						  -du[j+i*(M+1)])/(2.0*(j-1));
			}
			if(M%2==0){//Mo=Me+1 (extra odd mode)
				int j = 2*Mo-1;
				scro[(j-1)/2+i*Mo]=(du[j-1+i*(M+1)]
						    -du[j+1+i*(M+1)])/(2.0*j); 
			}
		}
		
		even_tri.solve(scre, nrhs);
		odd_tri.solve(scro, nrhs);

		/*
		 * unpack computed modes from scre and scro to du
		 */
		for(int i = 0; i < nrhs; i++){
			du[i*(M+1)] = 0;
			du[M+i*(M+1)] = 0;
			for(int j = 2; j <= 2*Me; j += 2){
				du[j+i*(M+1)] = scre[j/2-1+i*Me];
				du[j-1+i*(M+1)] = scro[j/2-1+i*Mo];
			}
			if(M%2==0){//Mo=Me+1 (extra odd mode)
				int j = 2*Mo-1;
				du[j+i*(M+1)] = scro[(j-1)/2+i*Mo];
			}
		}
	}
	/*
	 * from cheb series of du to that of u
	 */
	for(int i=0; i < nrhs; i++){
		u[i*(M+1)] = 0;
		u[M+i*(M+1)] = 0;
		for(int j=1; j < M; j++)
			u[j+i*(M+1)] = (du[j-1+i*(M+1)]
					- du[j+1+i*(M+1)])/(2.0*j);
	}

	/*
	 * transform u and du from cheb to space domain
	 */
	for(int i=0; i < nrhs; i++){
		trig->bwd(u+i*(M+1));
		trig->bwd(du+i*(M+1));
	}
}

/*
 * solvep() cut and paste and lightly modified
 */
void BVPSolvei::solvepx(const double *restrict f, double *restrict g,
			double *restrict u, double *restrict du, int nrhs){
	assrt(state==4);
	
	/*
	 * copy f to du
	 */
	for(int i=0; i < nrhs*(M+1); i++)
		du[i] = f[i];

	/*
	 * find cheb coeffs of f copied to du and of g
	 */
	for(int i=0; i < nrhs; i++){
		trig->fwd(du+i*(M+1));
		du[M+i*(M+1)] = 0;
		trig->fwd(g+i*(M+1));
		g[M+i*(M+1)] = 0;
	}
	
	{
		double *restrict scre = u;
		double *restrict scro = u + nrhs*Me;
		/*
		 * separate even odd modes of f in du
		 * integrate and pack into scre and scro
		 */
		for(int i=0; i < nrhs; i++){//Mo>=Me
			for(int j = 2; j <= 2*Me; j += 2){
				scre[j/2-1+i*Me]=(du[j-1+i*(M+1)]
						    -du[j+1+i*(M+1)])/(2.0*j)
					+ g[j+i*(M+1)];
				scro[j/2-1+i*Mo]=(du[j-2+i*(M+1)]
						  -du[j+i*(M+1)])/(2.0*(j-1))
					+ g[j-1+i*(M+1)];
			}
			if(M%2==0){//Mo=Me+1 (extra odd mode)
				int j = 2*Mo-1;
				scro[(j-1)/2+i*Mo]=(du[j-1+i*(M+1)]
						    -du[j+1+i*(M+1)])/(2.0*j)
					+ g[j+i*(M+1)]; 
			}
		}
		
		even_tri.solve(scre, nrhs);
		odd_tri.solve(scro, nrhs);

		/*
		 * unpack computed modes from scre and scro to du
		 */
		for(int i = 0; i < nrhs; i++){
			du[i*(M+1)] = 0;
			du[M+i*(M+1)] = 0;
			for(int j = 2; j <= 2*Me; j += 2){
				du[j+i*(M+1)] = scre[j/2-1+i*Me];
				du[j-1+i*(M+1)] = scro[j/2-1+i*Mo];
			}
			if(M%2==0){//Mo=Me+1 (extra odd mode)
				int j = 2*Mo-1;
				du[j+i*(M+1)] = scro[(j-1)/2+i*Mo];
			}
		}
	}
	/*
	 * from cheb series of du to that of u
	 */
	for(int i=0; i < nrhs; i++){
		u[i*(M+1)] = 0;
		u[M+i*(M+1)] = 0;
		for(int j=1; j < M; j++)
			u[j+i*(M+1)] = (du[j-1+i*(M+1)]
					- du[j+1+i*(M+1)])/(2.0*j);
	}

	/*
	 * transform u and du from cheb to space domain
	 */
	for(int i=0; i < nrhs; i++){
		trig->bwd(u+i*(M+1));
		trig->bwd(du+i*(M+1));
	}
}


