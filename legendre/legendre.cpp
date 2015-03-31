#include <cmath>
#include "../utils/utils.hh"
#include "../fft/trig.hh"
#include "legendre.hh"

int pow2(int n){
	while(n > 1){
		if (n%2 == 1)
			return 0;
		n  = n/2;
	}
	return 1;
}

FLTrans::FLTrans(int MM)
	:M(MM), dct(M)
{
	assrt(pow2(M) == 1 && M >= 64);
	wsave = new double[1l*199*375*M];
	plini_(&M, wsave);
}

FLTrans::~FLTrans(){
	delete[] wsave;
}

void FLTrans::fwd(double *v){
	dct.fwd(v);
	v[0] /= 2.0;
	v[M] = 0.0;
	plb_(&M, v, wsave);
	v[M] = 0.0;
}

void FLTrans::bwd(double *v){
	plf_(&M, v, wsave);
	v[0] *= 2.0;
	v[M] = 0.0;
	dct.bwd(v);
}
