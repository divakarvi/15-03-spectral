#ifndef __MKLTrig4July2013__
#define __MKLTrig4July2013__
#include <mkl.h>
#include <mkl_trig_transforms.h>
#include "../utils/utils.hh"

class DCT{//discrete cosine transform
private:
	int *ipar;
	double *dpar;
	int stat;
	DFTI_DESCRIPTOR_HANDLE handle;
public:
	DCT(int n);
	~DCT();
	void fwd(double *v){//size = n+1 doubles
		d_forward_trig_transform(v, &handle, ipar, dpar, &stat);
		assrt(stat==0);
	}
	void bwd(double *v){//size = n+1 doubles
		d_backward_trig_transform(v, &handle, ipar, dpar,&stat);
		assrt(stat==0);
	}
};

class DST{//discrete sine transform
private:
	int n;
	int *ipar;
	double *dpar;
	int stat;
	DFTI_DESCRIPTOR_HANDLE handle;
public:
	DST(int ni);
	~DST();
	void fwd(double *v){//size = n+1 doubles
		v[0] = v[n] = 0;
		d_forward_trig_transform(v, &handle, ipar, dpar, &stat);
		assrt(stat==0);
	}
	void bwd(double *v){//size = n+1 doubles
		v[0] = v[n] = 0;
		d_backward_trig_transform(v, &handle, ipar, dpar,&stat);
		assrt(stat==0);
	}
};

#endif
