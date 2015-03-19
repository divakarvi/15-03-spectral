#ifndef __dvutils25June2013__
#define __dvutils25June2013__
#include <cstdio>
#include <iostream>
#include <fstream>
const double PI = 3.1415926535897932384e+00;

/*
 * based on http://stackoverflow.com/questions/5252375/custom-c-assert-macro
 * works outside debug mode too
 */
#define assrt(cond) \  
    do \  
    { \  
	 if (!(cond)) \  
	 { \  
		  std::cout<<"ASSRT FAILED:"<<#cond<<std::endl	\
			   <<"file: "<<__FILE__<<std::endl	\
			   <<"line: "<<__LINE__<<std::endl;	\
		 __builtin_trap();				\
        } \  
    } while(0)  



 /*
  * v[i] = fabs(v[i]) for i=0...n-1
  */
void array_abs(double *v, int n);

/*
 * returns max of abs values of v[0..n-1]
 */
double array_max(double *v, int n);
 /*
  * display n entries of v on stdout
  */
void array_show(double *v, int n, const char* mesg);
/*
 * v = v - w
 */
void array_diff(double *restrict v, double *restrict w, int n);
 /*
  * copy v = w
  */
 void array_copy(double *restrict v, double *restrict w, int n);
/*
 * v[]   = array to be output (column major)
 * m     = num of rows
 * n     = num of cols 
 * fname = name of file for output
 */
void array_out(double *v, int m, int n, const char *fname);
 /*
  * same as above except for lda
  * i,j th entry = v[i+j*lda]
  */
void array_out(double *v, int m, int n, int lda, const char *fname);
/*
 * counterpart of array_out()
 * size = number of entries of v
 */
void array_in(double *v, int size,  const char* fname);

/*
 * A = 4 doubles in column major format
 * rhs = 2 doubles
 * x = soln of Ax=rhs
 */
void solve2x2(double *restrict A, double *restrict rhs, double *restrict x);


/*
 * reads x from fname
 */
void file2int(const char* fname, int& x);
void file2double(const char* fname, double& x);

/*
 * writes x to fname
 */
void int2file(int x, const char* fname);
void double2file(double x, const char* fname);

 /*
  * verify if the directory exists
  * if not create it
  */
void verify_dir(const char *dir);

/*
 * const char *dir = name of directory
 * const char *pfs = file name prefix
 * removes all files in directory starting with prefix
 */
void mop_dir(const char* dir, const char *pfx);

/*
 * prints the file to cout (omitting unprintable charaters)
 * the contents of the file are put inside a box
 */
void box_file(const char* fname, const char* mesg);

/*
 * getpid() and use it to copy /proc/pid/status to std::cout
 */
void print_proc_status(const char* mesg);
#endif
