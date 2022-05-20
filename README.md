# de_int_fi
This is the numerical integration program for Fourier-type integrals, using double exponential transformation. 
This code can calculate the sine type integral and the cosine type integral from a to infinity shown by the following formulae.
It also supports the integration of complex functions.
I made it to integrate a complex function which has a singular point.

$$ I = \int_a^{\infty} f(x) \sin(\omega x) \\,dx $$

$$ I = \int_a^{\infty} f(x) \cos(\omega x) \\,dx $$

$f(x)$ needs to be analytic within the integral region, and non-oscillating function.

$$ I = \int_a^{\infty} f(x) \\,dx, \\ \\ f(x) {\approx} g(x) \sin(\omega x), \\ \\ x{\gg}a $$

$$ I = \int_a^{\infty} f(x) \\,dx, \\ \\ f(x) {\approx} g(x) \cos(\omega x), \\ \\ x{\gg}a $$

$g(x)$ needs to be analytic within the integral region, and non-oscillating function.


## Usage of the example code  
 1. type 'make' command to compile.
 2. type './example.out' to run.

Please see de_int_fi.h for detail of functions, example.c for detail of function usages.


## Regarding the examples  
- Example 1  
  $$ I = \int_0^{\infty} \frac{\sin x}{x} \\,dx = \frac{\pi}{2} $$

- Example 2  
  $$ I = \int_0^{\infty} \frac{\cos x}{x^2+1} \\,dx = \frac{\pi}{2e^{}} $$
  
  $e^{}$ is Napier's constant ( Euler's number ).

- Example 3  
  $$ I = \int_0^{\infty} \frac{\sin x}{x+1+i} \\,dx = \mathrm{Ci}(z) \sin z + \left( \frac{\pi}{2} - \mathrm{Si}(z) \right) \cos z, \\ \\ z=1+i. $$

  $\mathrm{Ci}(z)$ is cosine integral function, $\mathrm{Si}(z)$ is sine integral function.
  
- Example 4  
  $$ I = \frac{2}{\pi i} \int_0^{\infty} \exp\left(-\frac{x-i}{\sqrt{2}}\right)\frac{\cos\left(\frac{\sqrt{x^2-2ix}}{\sqrt{2}}\right)}{\sqrt{x^2-2ix}} \\,dx = H_0^{(1)}(1) $$
  
  $H_0^{(1)}(z)$ is Hankel function 1st-kind order 0.


## Reference
Ooura, Takuya, and Masatake Mori. "A robust double exponential formula for Fourier-type integrals." Journal of computational and applied mathematics 112.1-2 (1999): 229-241.
