#ifndef __DVLEGENDRE31March2015__
#define __DVLEGENDRE31March2015__
#include "../fft/trig.hh"

/*
 * defined in legcheb.o
 */ 
extern "C"{
	void plini_(int *n, double *wsave);
	void plf_(int *n, double *v, double *wsave);
	void plb_(int *n, double *v, double *wsave);
};

/*
 * FLTrans = Fast Legendre Transformation
 */
class FLTrans{
private:
	int M;
	DCT dct;
	double *wsave;
public:
	FLTrans(int MM);
	~FLTrans();
	/*
	 *  input: v[0..M] of length M+1, v[i] = function value at cos(i*PI/M)
	 * output: v[0...M-1] are coeffs of Legendre polynomials
	 *         v[M] set to zero
	 */
	void fwd(double *v);
	/*
	 *  input: v[0..M] coeffs of Legendre polynomials, v[M] ignored
	 * output: v[0..M], v[i] = function value at cos(i*PI/M)
	 */
	void bwd(double *v);
	
};

#endif
