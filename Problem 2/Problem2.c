#include <stdio.h>
#include <math.h>
#include <conio.h>

// Definition:
// endpoint a, b, grid spacing h, initial condition alpha, n subintervals, exact_sol is exact solution
// Set endpoint a, b and initial condition alpha
// Following the given 0 <= t <= 2, y(0) = 0.5
const double a = 0, b = 2, alpha = 0.5;
double h, exact_sol = 0;

double dydt(double t, double y) {
  // dydt = y'(t, y) = y - t^2 + 1
  return y - pow(t, 2) + 1;
}

double y(double t) { 
  // 1 + t^2 + 2t - (1/2)*(e^t)
  return 1 + pow(t, 2) + 2*t - (0.5)*exp(t);
}

void plot_point(FILE *result, double a, double b, double step) {
  //Find exact solution
  while(a <= b) {
    fprintf(result, "%.10lf\t %.10lf\n", a, y(a));
    a += step;
  }
  exact_sol = y(a);
  printf("Exact solution: %.10lf\n", exact_sol);
}

double Absolute_error(double P, double p) {
  //P is the exact solution, p is approximate solution
  return fabs(P - p);
}

void Euler(FILE *result, double a, double b, double h, double alpha) {
  // Approximate solution using Euler method
  int n = (int)(b - a)/h;
  double t = a;
  double w = alpha;
  //printf("%.10lf\t %.10lf\n", t, w) //if you want to see all step;
  fprintf(result, "%.10lf\t %.10lf\n", t, w); 
  for(int i = 1; i <= n; i++) {
    w += h*dydt(t, w);
    t = a + i*h;
    //printf("%.10lf\t %.10lf\n", t, w);
    fprintf(result, "%.10lf\t %.10lf\n", t, w);
  }
  printf("Approximate solution: %.10lf\n", w);
  if(h == 0.01)
    printf("Absolute error: %.10lf\n", Absolute_error(exact_sol, w));
}

double deriv(double t, double y) {
  // first derivative of f(t, y(t)) = y(t) - t^2 + 1
  // with respect to the variable t
  // f'(t, y(t)) = y - t^2 + 1 -2t
  return y - pow(t, 2) + 1 - 2*t;
}

void Taylor2(FILE *result, double a, double b, double h, double alpha) {
  // Approximate solution using Taylor method of order 2
  int n = (int)(b - a)/h;
  double t = a;
  double w = alpha;
  //printf("%.10lf\t %.10lf\n", t, w);
  fprintf(result, "%.10lf\t %.10lf\n", t, w);
  for(int i = 1; i <= n; i++) {
    w += h*(dydt(t, w) + (h/2)*deriv(t, w));
    t = a + i*h;
    //printf("%.10lf\t %.10lf\n", t, w);
    fprintf(result, "%.10lf\t %.10lf\n", t, w);
  }
  printf("Approximate solution: %.10lf\n", w);
  if(h == 0.01)
    printf("Absolute error: %.10lf\n", Absolute_error(exact_sol, w));
}

void RK2(FILE *result, double a, double b, double h, double alpha) {
  // Approximate solution using Runge-Kutta method of order 2
  int n = (int)(b - a)/h;
  double K1, K2, t, w;
  t = a;
  w = alpha;
  //printf("%.10lf\t %.10lf\n", t, w);
  fprintf(result, "%.10lf\t %.10lf\n", t, w);
  for(int i = 1; i <= n; i++) {
    K1 = dydt(t, w);
    t = a + i*h;
    K2 = dydt(t, w + h*K1);
    w += (h/2)*(K1 + K2); 
    //printf("%.10lf\t %.10lf\n", t, w);
    fprintf(result, "%.10lf\t %.10lf\n", t, w);
  }
  printf("Approximate solution: %.10lf\n", w);
  if(h == 0.01)
    printf("Absolute error: %.10lf\n", Absolute_error(exact_sol, w));
}

void RK4(FILE *result, double a, double b, double h, double alpha) {
  // Approximate solution using Runge-Kutta method of order 4
  int n = (int)(b - a)/h;
  double t = a;
  double w = alpha;
  //printf("%.10f\t %.10lf\n", t, w);
  fprintf(result, "%.10lf\t %.10lf\n", t, w);
  for(int i = 1; i <= n; i++) {
    double K1 = h*dydt(t, w);
    double K2 = h*dydt(t + h/2, w + K1/2);
    double K3 = h*dydt(t + h/2, w + K2/2);
    double K4 = h*dydt(t + h, w + K3);
    w += (K1 + 2*K2 + 2*K3 + K4)/6;
    t = a + i*h;
    //printf("%.10lf\t %.10lf\n", t, w);
    fprintf(result, "%.10lf\t %.10lf\n", t, w);
  }
  printf("Approximate solution: %.10lf\n", w);
  if(h == 0.01)
    printf("Absolute error: %.10lf\n", Absolute_error(exact_sol, w));
}

int main(void) {
  FILE *fp = NULL;
  fp = freopen("Exact solution.txt", "w", stdin);
  plot_point(fp, 0, 2, 0.01);
  //Exact solution 5.3054719505
  
  puts("");
  h = 0.1;
  puts("Compute approximate solution using numerical method");
  printf("with h = 0.1\n");

  fp = freopen("Euler h01.txt", "w", stdin);
  printf("Euler method \n");
  Euler(fp, a, b, h, alpha);
  
  puts("");
  fp = freopen("Taylor2 h01.txt", "w", stdin);
  printf("Taylor of order 2 \n");
  Taylor2(fp, a, b, h, alpha);
  
  puts("");
  fp = freopen("RungeKutta4_h01.txt", "w", stdin);
  printf("Runge-Kutta of order 4\n");
  RK4(fp, a, b, h, alpha);
  
  puts("");
  fp = freopen("RungeKutta2 h01.txt", "w", stdin);
  printf("Runge-Kutta of order 2\n");
  RK2(fp, a, b, h, alpha);

  puts("");
  puts("==================================");
  h = 0.01;
  puts("Compute with h = 0.01");

  fp = freopen("Euler h001.txt", "w", stdin);
  printf("Euler method \n");
  Euler(fp, a, b, h, alpha);
  
  puts("");
  fp = freopen("Taylor2 h001.txt", "w", stdin);
  printf("Taylor of order 2 \n");
  Taylor2(fp, a, b, h, alpha);
  
  puts("");
  fp = freopen("RungeKutta4 h001.txt", "w", stdin);
  printf("Runge-Kutta of order 4 \n");
  RK4(fp, a, b, h, alpha);
  
  puts("");
  fp = freopen("RungeKutta2 h001.txt", "w", stdin);
  printf("Runge-Kutta of order 2 \n");
  RK2(fp, a, b, h, alpha);
  
  fclose(fp);
  getch();
} 
