// numerical integration program for fourier type integral using double exponential formula.
#ifndef DE_INT_FI_H_
#define DE_INT_FI_H_

#include <math.h>
#include <complex.h>

#define eerc 5.0 // parameter for error estimation
#define LMAX 256 // limit of summation

double deintd_fts(double (*f)(double,void *),void *pa,double a,double omega,double eps,int *erc);
double deintd_ftc(double (*f)(double,void *),void *pa,double a,double omega,double eps,int *erc);

double deintd_ffs(double (*f)(double,void *),void *pa,double a,double omega,double eps,int *erc);
double deintd_ffc(double (*f)(double,void *),void *pa,double a,double omega,double eps,int *erc);
/*  
   deint_fts : integrate over (a,infinity), integrand f(x) * sin(omega x), f(x) is not a oscillation function.
   deint_ftc : integrate over (a,infinity), integrand f(x) * cos(omega x), f(x) is not a oscillation function.
   
   deint_ffs : integrate over (a,infinity), integrand f(x) include oscillation part. 
               x>>a f(x) ~ g(x) sin(omega x), g(x) is not a oscillation function.
   deint_ffc : integrate over (a,infinity), integrand f(x) include oscillation part. 
               x>>a f(x) ~ g(x) cos(omega x), g(x) is not a oscillation function.

   -- arguments -- 
   double (*f)(double,void *) : function pointer of integrand
                                f(variable of integration, pointer of paramter list casted to void type)
   void *pa                   : pointer of parameter list 
   double a                   : lower limit of integration
   double omega               : angular frequency
   double eps                 : relative error requested
   int *erc                   : error code 
                                erc = 0 : normal termination
                                erc =-1 : upward integration terminated by error
                                erc =-2 : downward integration terminated by error
*/

// for a integrand which returns a complex value
double complex deintz_ffs(double complex (*f)(double,void *),void *pa,double a,double omega,double eps,int *erc);
double complex deintz_ffc(double complex (*f)(double,void *),void *pa,double a,double omega,double eps,int *erc);

#endif /* DE_INT_FI_H_ */
