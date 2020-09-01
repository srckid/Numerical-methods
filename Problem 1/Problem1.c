#include <stdio.h>
#include <math.h>
#include <conio.h>

double f(double x) {
  // Declare f(x) = e^x - 4(x^3) + 1
  return exp(x) - 4*(x*x*x) + 1;
}

void func_point(FILE *fp, double start, double end, double step) {
  // Save function point x, y to the text file
  // use (x,y) to plot the graph
  double x, y;
  for(x = start; x <= end; x += step) {
    y = f(x);
    fprintf(fp, "%lf\t %lf\n", x, y);
  }
}

// Definition: initial p0, endpoint a, b, tolerance TOL, n iteration 

// Find root using Bisection method
void Bisection(double a, double b, double TOL, int n) {
  if(f(a)*f(b) <= 0) {
    int i = 1;
    double c;
    do {
      // Set midpoint c
      c = (a + b)/2;
      if(f(a)*f(c) > 0) {
        a = c;
      }
      else if(f(a)*f(c) < 0) {
        b = c;
      }
      i++;     
    } while(fabs(a - b) >= TOL && i <= n);
    // Stop when |a - b| < tolerance and i > n 
    printf("%.8lf", c);
    return;
  }
  else
    printf("The method failed after %d iteration", n);
}

double g(double x) {
  // Declare g(x) = ((e^x + 1)/4)^(1/3)
  // cdrt is the function to return cube root of argument 
  return cbrt((exp(x)+1)/4);
}
// Find root using Fixed Point iteration method
void Fixed_Point(double p0, double TOL, int n) {  
  int i = 1;
  double p;
  while(i <= n) {
    p = g(p0);
    if(fabs(p - p0) < TOL) {
      printf("%.8lf", p);
      return;
    }
    i++;
    p0 = p;
  }
  printf("The method failed after %d iteration", n);
}

double F(double x) {
  // Declare the derivation F(x) = f'(x) = e^x - 12(e^2) 
  return exp(x) - 12*pow(x,2);
}
// Find root using Newton Rapshson method
void Newton(double p0, double TOL, int n) { 
  int i = 1;
  while(i <= n) {
    double p = p0 - f(p0)/F(p0);
    if(fabs(p - p0) < TOL) {
      printf("%.8lf\n", p);
      return;
    }
    i += 1;
    p0 = p;
  }
  printf("The method failed after %d iteration", n);
}

int main(void) {
  FILE *fp = NULL;
  // Make graph point and save into text file 
  // with interval[-5, 5] and step 0.1 to plot to find approximate solution
  fp = fopen("function_point.txt", "w");
  func_point(fp, -5, 5, 0.1);
  fclose(fp);
  // Estimate root interval through graph [a, b]
  // By graph, we estimate the root interval a = -2, b = 2 
  const double a = -2, b = 2, error = 0.0001, n = 10;
  // Find p10 with initial guess = 1, error, 10 iteration
  printf("Approximate p10\n");
  printf("using Bisection method     : "); Bisection(a, b, error, n);
  puts("");
  printf("using Fixed Point iteration: "); Fixed_Point(1, error, n);
  puts("");
  printf("using Newton-Raphson method: "); Newton(1, error, n);
  
  getch();
  return 0;
}
