#include <stdio.h>
#include <math.h>
#include <conio.h>

#define A 0
#define B 1
#define H 0.00001

// Definition:
// endpoint a b, grid spacing h, n subintervals

double f(double x) {
  // Declare function f(x) = 1/(1+x*x)
  return 1 / (1 + x * x);
}

double Composite_Trapezoidal(double a, double b, double h) {
  // Approximate value using Composite Tranpezoidal rule
  int n = (b - a) / h;
  // Computing sum of first and last terms  
  double s = f(a) + f(b); 
  // Adding middle terms in above formula 
  for (int i = 1; i < n; i++) 
    s += 2 * f(a + i*h); 
  // h/2 indicates (b-a)/2n. Multiplying h/2 with s. 
  return (h / 2) * s; 
}

double Composite_Simpson(double a, double b, double h) {
  // Approximate value using Composite Simpson's rule 
  int n = (b - a) / h;
  // Set summation
  // s0 is summation of first and last terms
  // s1 is summation of f(X(2i-1))
  // s2 is summation of f(X(2i))
  double s0 = f(a) + f(b);
  double s1 = 0;
  double s2 = 0;
  for(int i = 1; i < n; i++) {
    double X = a + i*h;
    // if i is even then set s2 = s2 + f(x)
    // else set s1 = s1 + f(x)
    if(i % 2 == 0)
      s2 += f(X);
    else
      s1 += f(X);
  }
  // Compute final result
  return (h / 3) * (s0 + 2*s2 + 4*s1);
}

// Approximate pi number through pi/4 pass by integral [0,1] of f(x) = 1/(1+x*x)
#define pi(result) result * 4

int main(void) {
  printf("Composite Trapezoidal: %.8lf", Composite_Trapezoidal(A, B, H));
  puts("");
  printf("Composite_Simpson: %.8lf", Composite_Simpson(A, B, H));
  puts("");
  printf("Approximate pi: %.8lf %.8lf", pi(Composite_Simpson(A, B, H)), pi(Composite_Trapezoidal(A, B, H)));
  
  getch();
  return 0;
}
