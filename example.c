#include <stdio.h>
#include "de_int_fi.h"

// parameter set
typedef struct param{
  double omega; // variable for pass the angular frequency to the integrand
  int cc;       // counter of function calls

}PM;

// function prototype of the integrand
// Example 1, sin type function
double s1(double x,void *pa);         // integrand for deintd_fts ( oscillation part is excluded )
double s2(double x,void *pa);         // integrand for deintd_ffs ( oscillation part is included )
// Example 2, cos type function 
double c1(double x,void *pa);         // integrand for deintd_ftc ( oscillation part is excluded )
double c2(double x,void *pa);         // integrand for deintd_ffc ( oscillation part is included )
// Example 3, sin type function
double complex zs(double x,void *pa); 
// Example 4, cos type function 
double complex zc(double x,void *pa); 

int main()
{
  PM pa;
  double complex Iz;
  double Id;
  double a=0.0;
  double omega=1.0;
  double eps=1.0e-15;
  int erc;

  pa.omega=omega;
  printf("-- Example 1 ( sin type function ) --\n");
  printf("I    = int_0^inf sin(x)/x dx\n");
  printf("     = %15.14e ( = pi/2 )\n",M_PI/2.0);
  pa.cc=0;
  Id=deintd_fts(s1,&pa,a,omega,eps,&erc);
  printf("I_fts= %15.14e, erc = %d, function call %d\n",Id,erc,pa.cc);
  pa.cc=0;
  Id=deintd_ffs(s2,&pa,a,omega,eps,&erc);
  printf("I_ffs= %15.14e, erc = %d, function call %d\n\n",Id,erc,pa.cc);

  printf("-- Example 2 ( cos type function ) --\n");
  printf("I    = int_0^inf cos(x)/(x^2+1) dx\n");
  printf("     = %15.14e ( = pi/(2*e) )\n",M_PI/(2.0*exp(1))); 
  pa.cc=0;
  Id=deintd_ftc(c1,&pa,a,omega,eps,&erc);
  printf("I_ftc= %15.14e, erc = %d, function call %d\n",Id,erc,pa.cc);
  pa.cc=0;
  Id=deintd_ffc(c2,&pa,a,omega,eps,&erc);
  printf("I_ffc= %15.14e, erc = %d, function call %d\n\n",Id,erc,pa.cc);

  printf("-- Example 3 --\n");
  printf("I    = int_0^inf sin(x)/(x+1+i) dx\n");
  printf("     = Ci(z)*sin(z)+(pi/2-Si(z))*cos(z), z=1+i. Ci(z) : cosine integral function, Si(z) : sine integral function\n");
  printf("     = 4.79410254818214e-01 -2.63993546979857e-01 i \n");
  pa.cc=0;
  Iz=deintz_ffs(zs,&pa,a,1.0,eps,&erc);
  printf("I_ffs= %15.14e %+15.14e i, erc = %d, function call %d\n\n",creal(Iz),cimag(Iz),erc,pa.cc);

  printf("-- Example 4 --\n");
  printf("I    = 2/(pi*i) * int_0^inf exp(-(x-i)/sqrt(2))*cos(sqrt(x^2-2*i*x)/sqrt(2))/sqrt(x^2-2*i*x) dx\n");
  printf("     = H_0^(1)(1), H_0^(1)(x) : hankel function 1st-kind order 0\n");
  printf("     = 7.65197686557966e-01 +8.82569642156769e-02 i\n");
  pa.cc=0;
  Iz=2.0/(M_PI*I)*deintz_ffc(zc,&pa,a,1.0/sqrt(2.0),eps,&erc);
  printf("I_ffc= %15.14e %+15.14e i, erc = %d, function call %d\n",creal(Iz),cimag(Iz),erc,pa.cc);

  return 0;
}

// Example 1, oscillation part is excluded.
double s1(double x,void *pa)
{
  PM *t=(PM*)pa;
  t->cc+=1;
  return 1.0/x;
}

// Example 1, oscillation part is included.
double s2(double x,void *pa)
{
  PM *t=(PM*)pa;
  t->cc+=1;
  return 1.0/x*sin(t->omega*x);
}

// Example 2, oscillation part is excluded.
double c1(double x,void *pa)
{
  PM *t=(PM*)pa;
  t->cc+=1;
  return 1.0/(x*x+1.0);
}

// Example 2, oscillation part is included.
double c2(double x,void *pa)
{
  PM *t=(PM*)pa;
  t->cc+=1;
  return 1.0/(x*x+1.0)*cos(t->omega*x);
}

// Example 3
double complex zs(double x,void *pa)
{
  PM *t=(PM*)pa;
  t->cc+=1;
  return sin(x)/(x+1.0+1.0*I);
}

// Example 4
double complex zc(double x,void *pa)
{
  PM *t=(PM*)pa;
  double complex sr,cc;
  t->cc+=1;

  sr=csqrt(x*x-2.0*I*x);
  cc=ccos(sr/sqrt(2.0));
  return cexp(-(x-I)/sqrt(2.0))*cc/sr;
}

