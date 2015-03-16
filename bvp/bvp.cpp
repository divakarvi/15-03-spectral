#include "bvp.hh"
#include "../utils/utils.hh"

BVPSolve::BVPSolve(double aa, int MMi, int MMo)
	:state(0), a(aa), Mi(MMi), Mo(MMo),
	 banded(2*Mo)
{
	assrt(Mo >= 1); 
	assrt(state == 0);
	state = 1;

	typedef BVPSolvei *bvpsolveptr;
	bvpi = new bvpsolveptr[Mo];
}

BVPSolve::~BVPSolve(){
	assrt(state == 5);

	for(int i=0; i < Mo; i++)
		delete bvpi[i];
	delete[] bvpi;
}

void BVPSolve::init_oy(const double *oyi){
	assrt(state == 1);
	state = 2;
	oy = oyi;

	for(int i=0; i < Mo; i++){
		double s = 2.0/(oy[i] - oy[i+1]);
		bvpi[i] = new BVPSolvei(a/s, Mi);
	}
}

void BVPSolve::initcheb(double *chebi){
	assrt(state==2);
	state=3;
	for(int i=0; i < Mo; i++)
		bvpi[i]->initcheb(chebi);
}

void BVPSolve::initdct(DCT& dcti){
	assrt(state==3);
	state = 4;
	for(int i=0; i < Mo; i++){
		bvpi[i]->initdct(dcti);
		bvpi[i]->inithmg();
	}
}

void BVPSolve::initbanded(){
	assrt(state==4);
	state = 5;
	if(Mo > 1){
		double tmp[2*Mo+4];
		for(int i=0; i < Mo-1; i++){
			double h2p1 = bvpi[i+1]->geth2()[0];
			tmp[2*i] = 0;
			tmp[2*i+1] = -h2p1;
		}
		banded.setu2(tmp);

		double h2p1_0 = bvpi[0]->geth2()[0];
		tmp[0] = -h2p1_0;
		for(int i=0; i < Mo-1; i++){
			double h1p1 = bvpi[i+1]->geth1()[0];
			double d2p1 = bvpi[i+1]->getd2()[0];
			double scale_p = 2.0/(oy[i+1] - oy[i+2]);
			tmp[2*i+1] = -h1p1;
			tmp[2*i+2] = -d2p1 * scale_p;
		}
		banded.setu1(tmp);

		double h1p1_0 = bvpi[0]->geth1()[0];
		tmp[0] = -h1p1_0;
		for(int i=0; i < Mo; i++){
			double h2m1 = bvpi[i]->geth2()[Mi];
			double d1p1;
			if(i+1 < Mo)
				d1p1 = bvpi[i+1]->getd1()[0];
			double scale_p = 2.0/(oy[i+1] - oy[i+2]);
			tmp[2*i+1] = h2m1;
			tmp[2*i+2] = -d1p1 * scale_p;
		}
		banded.setdiag(tmp);

		for(int i=0; i < Mo; i++){
			double h1m1 = bvpi[i]->geth1()[Mi];
			double d2m1 = bvpi[i]->getd2()[Mi];
			double scale_m = 2.0/(oy[i] - oy[i+1]);
			tmp[2*i] = h1m1;
			tmp[2*i+1] = d2m1 * scale_m;
		}
		banded.setl1(tmp);

		for(int i=0; i < Mo-1; i++){
			double d1m1 = bvpi[i]->getd1()[Mi];
			double scale_m = 2.0/(oy[i] - oy[i+1]);			
			tmp[2*i] = d1m1 * scale_m;
			tmp[2*i+1] = 0;
		}
		banded.setl2(tmp);
		banded.factor();
	}
	else{
		double h1p1 = bvpi[0]->geth1()[0];
		double h2p1 = bvpi[0]->geth2()[0];
		double h1m1 = bvpi[0]->geth1()[Mi];
		double h2m1 = bvpi[0]->geth2()[Mi];
		B[0] = -h1p1; B[2] = -h2p1;
		B[1] =  h1m1; B[3] =  h2m1;
	}
		
}

void BVPSolve::solve(const double *f, 
		     double *restrict u, double *du,
		     double plusone, double minusone,
		     double *restrict g){
	assrt(state==5);
	
	double tmpf[(Mi+1)*Mo];
	for(int i=0; i < Mo; i++){
		double scale = 2.0/(oy[i] - oy[i+1]);
		/*
		 * scale to transform domain
		 */
		for(int j=0; j <= Mi; j++) 
			tmpf[j + i*(Mi+1)] = f[j + i*(Mi+1)]/scale/scale;

		if (g != NULL){
			for(int j=0; j <= Mi; j++)
				g[j+i*(Mi+1)] /= scale;
		}

		/*
		 * solve bvp
		 */
		if(g != NULL)
			bvpi[i]->solvepx(tmpf + i*(Mi+1), 
					 g + i*(Mi+1),
					 u + i*(Mi+1),
					 du + i*(Mi+1),
					 1);
		else
			bvpi[i]->solvep(tmpf + i*(Mi+1), 
					u + i*(Mi+1),
					du + i*(Mi+1),
					1);
		/*
		 * scale derv back to original domain
		 */
		for(int j=0; j <= Mi; j++)
			du[j + i*(Mi+1)] *= scale;
	}
	
	int M = Mi;
	double rhs[2*Mo+2];
	rhs[0] = -plusone+u[0];
	rhs[2*Mo-1] = minusone - u[M+(Mo-1)*(M+1)];
	for(int m=1; m < Mo; m++){
		rhs[2*m-1] = u[m*(M+1)] - u[M+(m-1)*(M+1)];
		rhs[2*m] = du[m*(M+1)] - du[M+(m-1)*(M+1)];
	}

	if(Mo > 1)
		banded.solve(rhs, 1); //nrhs = 1
	else
		solve2x2(B, rhs, rhs+2);

	double *coeffs = (Mo>1)?rhs:rhs+2;
	for(int i=0; i < Mo; i++){
		double *h1 = bvpi[i]->geth1();
		double *h2 = bvpi[i]->geth2();
		double *d1 = bvpi[i]->getd1();
		double *d2 = bvpi[i]->getd2();
		double scale = 2.0/(oy[i] - oy[i+1]);
		for(int j=0; j <= Mi; j++){
			u[j+i*(Mi+1)] += coeffs[2*i]*h1[j] 
				+ coeffs[2*i+1]*h2[j];
			du[j+i*(Mi+1)] += coeffs[2*i] * d1[j] * scale 
				+ coeffs[2*i+1] * d2[j] * scale;

		}
	}
}

void BVPSolve::hmg10n01(double *restrict pa, double *restrict dpa,
			   double *restrict pb, double *restrict dpb){
	assrt(state==5);
	int nrhs = Mo;
	
	double plusone = 0;
	double minusone = 1;

	int M = Mi;
	double rhs[2*Mo+2];
	rhs[0] = -plusone;
	rhs[2*Mo-1] = minusone;
	for(int m=1; m < Mo; m++){
		rhs[2*m-1] = 0;
		rhs[2*m] = 0;
	}

	if(Mo > 1)
		banded.solve(rhs, 1); //nrhs = 1
	else
		solve2x2(B, rhs, rhs+2);
	
	double *coeffs = (Mo>1)?rhs:rhs+2;
	for(int i=0; i < Mo; i++){
		double *h1 = bvpi[i]->geth1();
		double *h2 = bvpi[i]->geth2();
		double *d1 = bvpi[i]->getd1();
		double *d2 = bvpi[i]->getd2();
		double scale = 2.0/(oy[i] - oy[i+1]);
		for(int j=0; j <= M; j++){
			pb[j+i*(M+1)] = coeffs[2*i]*h1[j] 
				+ coeffs[2*i+1]*h2[j];
			dpb[j+i*(M+1)] = coeffs[2*i] * d1[j] * scale
				+ coeffs[2*i+1] * d2[j] * scale;
		}
	}
	
	plusone = 1;
	minusone = 0;
	rhs[0] = -plusone;
	rhs[2*Mo-1] = minusone;
	for(int m=1; m < Mo; m++){
		rhs[2*m-1] = 0;
		rhs[2*m] = 0;
	}

	if(Mo > 1)
		banded.solve(rhs, 1); //nrhs = 1
	else
		solve2x2(B, rhs, rhs+2);
	coeffs = (Mo>1)?rhs:rhs+2;
	for(int i=0; i < Mo; i++){
		double *h1 = bvpi[i]->geth1();
		double *h2 = bvpi[i]->geth2();
		double *d1 = bvpi[i]->getd1();
		double *d2 = bvpi[i]->getd2();
		double scale = 2.0/(oy[i] - oy[i+1]);
		for(int j=0; j <= M; j++){
			pa[j+i*(M+1)] = coeffs[2*i] * h1[j] 
				+ coeffs[2*i+1] * h2[j];
			dpa[j+i*(M+1)] = coeffs[2*i] * d1[j] * scale
				+ coeffs[2*i+1] * d2[j] * scale;
		}
	}
}
