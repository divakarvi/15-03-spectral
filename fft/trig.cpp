#include "trig.hh"
#include "../utils/utils.hh"

DCT::DCT(int n){
	int tt_type;
	tt_type = MKL_COSINE_TRANSFORM;
	ipar = new int[128];
	dpar = new double[5*n/2+2];
	d_init_trig_transform(&n, &tt_type, ipar, dpar, &stat);
	double *f = (double *)MKL_malloc((n+1)*sizeof(double), 32);
	d_commit_trig_transform(f, &handle, ipar, dpar, &stat);
	MKL_free(f);
}

DCT::~DCT(){
	free_trig_transform(&handle, ipar, &stat);
	delete[] ipar;
	delete[] dpar;
}

DST::DST(int ni)
	:n(ni){
	int tt_type;
	tt_type = MKL_SINE_TRANSFORM;
	ipar = new int[128];
	dpar = new double[5*n/2+2];
	d_init_trig_transform(&n, &tt_type, ipar, dpar, &stat);
	double *f = (double *)MKL_malloc((n+1)*sizeof(double), 32);
	f[0] = f[n] = 0;
	d_commit_trig_transform(f, &handle, ipar, dpar, &stat);
	MKL_free(f);
}

DST::~DST(){
	free_trig_transform(&handle, ipar, &stat);
	delete[] ipar;
	delete[] dpar;
}
