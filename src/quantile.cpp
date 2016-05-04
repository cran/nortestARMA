// g++ -c calcstat.cpp -o calcstat.o -I"/usr/lib/R/include"
// g++ -shared -o calcstatC.so calcstat.o -I"/usr/lib/R/include" -L"/usr/lib" -lR

#include <R.h>
#include "Rmath.h"
#include <R_ext/Applic.h>

extern "C" {



  void pvalst(double* x, int *xlen, int* n, int* d, int* K, double *res, double* term1, double* term2, double* term3, double* term4) {
    
    double Delta(double x, int d);
    double delta(double x, int d);

    void f1t2(double *z, int n, void *ex2); // f1v1

    void f2t1(double *z, int n, void *ex1); // f1v2
    void f2t2(double *z, int n, void *ex2); // f2v2
    void f2t3af3t3a(double *z, int n, void *ex3); // f4v2
    void f2t3bf3t3b(double *z, int n, void *ex3); // f5v2

    void f3t1inside(double *z, int n, void *ex2); // f1v3
    void f3t1(double *y, int n, void *ex1); // f2v3
    void f3t2a(double *z, int n, void *ex3); // f3v3
    void f3t2b(double *z, int n, void *ex1); // f4v3

    void f3t4a(double *x1, int n, void *ex3);
    void f3t4b(double *x1, int n, void *ex3);
    void f3t4c(double *x1, int n, void *ex3);
    void f3t4d(double *x1, int n, void *ex3);
    void f3t4e(double *x1, int n, void *ex3);

    void Rdqags(integr_fn f, void *ex, double *a, double *b,
		double *epsabs, double *epsrel,
		double *result, double *abserr, int *neval, int *ier,
		int *limit, int *lenw, int *last,
		int *iwork, double *work);
    
    
    double terme1, terme2, terme3, terme4, somme;
    
    double *a, *b, *epsabs, *epsrel, *result, *ex1, *ex2, *ex3, *abserr, *work;
    int *last, *limit, *lenw, *ier, *neval, *iwork, i;
    
    epsabs = new double[1];
    epsrel = new double[1];
    ex1 = new double[1];
    ex2 = new double[2];
    ex3 = new double[3];
    limit = new int[1];
    lenw = new int[1];
    a = new double[1];
    b = new double[1];
    result = new double[1];
    abserr = new double[1];
    last = new int[1];
    ier = new int[1];
    neval = new int[1];
    

    ex2[1] = (double)d[0];
    ex3[0] = (double)n[0];
    ex3[2] = (double)d[0];
    epsabs[0] = 0.0001220703125;
    epsrel[0] = 0.0001220703125;


    limit[0] = 100;
    lenw[0] = 4 * limit[0];
    somme = 0.0;

    iwork = new int[limit[0]];
    work = new double[lenw[0]];

        
      
    if (K[0] == d[0] - 1) { //s=0: Une chi2(d) // Mais c'est impossible ...
      terme1 = 0.0;
      terme2 = 0.0;
      terme3 = 0.0;
      terme4 = 0.0;
      for (i = 1; i <= xlen[0]; i++) {
	somme = Delta(x[i - 1], d[0]);
	term1[0] = terme1;
	term2[0] = terme2;
	term3[0] = terme3;
	term4[0] = terme4;    
	res[i - 1] = 1.0 - somme;
      }
    }
      
    if (K[0] == d[0]) { // s=1 : Moi avec Ducharme
      terme3 = 0.0;
      terme4 = 0.0;
      for (i = 1; i <= xlen[0]; i++) {

	terme1 = Delta(x[i - 1], d[0]) * Delta(log(n[0]), 1);
	
	if (x[i - 1] <= log(n[0])) terme2 = 0.0;
	else { 	  
	  ex2[0] = x[i - 1];
	  a[0] = log(n[0]);
	  b[0] = x[i - 1];	  
	  Rdqags(f1t2, ex2, a, b, epsabs, epsrel, result, abserr, neval, ier, limit, lenw, last, iwork, work); // f1v1
	  terme2 = result[0];
	}
	
	somme = terme1 + terme2;
	term1[0] = terme1;
	term2[0] = terme2;
	term3[0] = terme3;
	term4[0] = terme4;    
	res[i - 1] = 1.0 - somme;
      }
    }
      
    if (K[0] == d[0] + 1) { // s=2 : Ducharme avec Fontez (voir A Smooth Test of Goodness-of-Fit for Growth Curves and Monotonic Nonlinear Regression Models, 2004
      terme4 = 0.0;
      for (i = 1; i <= xlen[0]; i++) {

	ex1[0] = (double)n[0];
	a[0] = 0.0;
	b[0] = log(n[0]);
	
	Rdqags(f2t1, ex1, a, b, epsabs, epsrel, result, abserr, neval, ier, limit, lenw, last, iwork, work);    // f1v2
	terme1 = Delta(x[i - 1], d[0]) * result[0]; // Son premier terme (p. 5)
	
	if (x[i - 1] <= log(n[0])) terme2 = 0.0;
	else { 
	  ex2[0] = x[i - 1];
	  a[0] = 0.0;
	  b[0] = x[i - 1] - log(n[0]);
	  Rdqags(f2t2, ex2, a, b, epsabs, epsrel, result, abserr, neval, ier, limit, lenw, last, iwork, work); // f2v2
	  terme2 = Delta(log(n[0]), 1) * (result[0] - Delta(log(n[0]), 1) * Delta(x[i - 1] - log(n[0]), d[0])); // Son 3eme terme (p. 5)
	}
	
	if (x[i - 1] <= 2 * log(n[0])) terme3 = 0.0;
	else { 
	  ex3[1] = x[i - 1];
	  a[0] = 0.0;
	  b[0] = log(n[0]);
	  Rdqags(f2t3af3t3a, ex3, a, b, epsabs, epsrel, result, abserr, neval, ier, limit, lenw, last, iwork, work); // Son 2eme terme (p. 5)  // f4v2
	  terme3 = result[0]; 
	  a[0] = log(n[0]);
	  b[0] = x[i - 1] - log(n[0]);
	  Rdqags(f2t3bf3t3b, ex3, a, b, epsabs, epsrel, result, abserr, neval, ier, limit, lenw, last, iwork, work); // Son 4eme terme (p. 5)  // f5v2
	  terme3 = terme3 + result[0];
	}
	
	somme = terme1 + terme2 + terme3;
	term1[0] = terme1;
	term2[0] = terme2;
	term3[0] = terme3;
	term4[0] = terme4;    
	res[i - 1] = 1.0 - somme;

      }
	
    }
    
    if (K[0] > d[0] + 1) { // s=3 : Moi tout seul
      for (i = 1; i <= xlen[0]; i++) {

	ex1[0] = (double)n[0];
	a[0] = 0;
	b[0] = log(n[0]);      
	Rdqags(f3t1, ex1, a, b, epsabs, epsrel, result, abserr, neval, ier, limit, lenw, last, iwork, work); // f2v3
	terme1 = Delta(x[i - 1], d[0]) * result[0];

	if (x[i - 1] <= log(n[0])) terme2 = 0.0;
	else { 
	  ex3[1] = x[i - 1];
	  a[0] = log(n[0]);
	  b[0] = x[i - 1];
	  Rdqags(f3t2a, ex3, a, b, epsabs, epsrel, result, abserr, neval, ier, limit, lenw, last, iwork, work); // f3v3
	  terme2 = result[0];

	  ex1[0] = (double)n[0];
	  a[0] = 0;
	  b[0] = log(n[0]);
	  Rdqags(f3t2b, ex1, a, b, epsabs, epsrel, result, abserr, neval, ier, limit, lenw, last, iwork, work); // f4v3
	  terme2 = terme2 * result[0];

	}
      
	if (x[i - 1] <= 2.0 * log(n[0])) terme3 = 0.0;
	else { 
	  ex3[1] = x[i - 1];
	  a[0] = 0;
	  b[0] = log(n[0]);
	  Rdqags(f2t3af3t3a, ex3, a, b, epsabs, epsrel, result, abserr, neval, ier, limit, lenw, last, iwork, work); // f4v2
	  terme3 = result[0];
	  a[0] = log(n[0]);
	  b[0] = x[i - 1] - log(n[0]);
	  Rdqags(f2t3bf3t3b, ex3, a, b, epsabs, epsrel, result, abserr, neval, ier, limit, lenw, last, iwork, work); // f5v2
	  terme3 = terme3 + result[0];
	  terme3 = Delta(log(n[0]), 1) * terme3;
	}

	if (x[i - 1] <= 3.0 * log(n[0])) terme4 = 0.0;
	else { 
	  ex3[1] = x[i - 1];
	  a[0] = 0.0;
	  b[0] = log(n[0]);
	  Rdqags(f3t4a, ex3, a, b, epsabs, epsrel, result, abserr, neval, ier, limit, lenw, last, iwork, work);
	  terme4 = result[0];
	  a[0] = 0.0;
	  b[0] = log(n[0]);
	  Rdqags(f3t4b, ex3, a, b, epsabs, epsrel, result, abserr, neval, ier, limit, lenw, last, iwork, work);
	  terme4 = terme4 + result[0];
	  a[0] = 0.0;
	  b[0] = log(n[0]);
	  Rdqags(f3t4c, ex3, a, b, epsabs, epsrel, result, abserr, neval, ier, limit, lenw, last, iwork, work);
	  terme4 = terme4 + result[0];
	  a[0] = log(n[0]);
	  b[0] = x[i - 1] - 2.0 * log(n[0]);
	  Rdqags(f3t4d, ex3, a, b, epsabs, epsrel, result, abserr, neval, ier, limit, lenw, last, iwork, work);
	  terme4 = terme4 + result[0];
	  a[0] = log(n[0]);
	  b[0] = x[i - 1] - 2.0 * log(n[0]);
	  Rdqags(f3t4e, ex3, a, b, epsabs, epsrel, result, abserr, neval, ier, limit, lenw, last, iwork, work);
	  terme4 = terme4 + result[0];
	}

	somme = terme1 + terme2 + terme3 + terme4;
	term1[0] = terme1;
	term2[0] = terme2;
	term3[0] = terme3;
	term4[0] = terme4;    
	res[i - 1] = 1.0 - somme;

      }
    }
      
      
  

  //On libere de la memoire
  delete[] epsabs;
  delete[] epsrel;
  delete[] ex1;
  delete[] ex2;
  delete[] ex3;
  delete[] limit;
  delete[] lenw;
  delete[] a;
  delete[] b;
  delete[] result;
  delete[] abserr;
  delete[] last;
  delete[] ier;
  delete[] neval;
  delete[] iwork;
  delete[] work;


}
      
      
  

  double pnorm(double x, double mu, double sigma, int lower_tail,int give_log);
  void R_cumsum(double *x, int *n, double *na_value, double *ans);
  
  typedef void integr_fn(double *x, int n, void *ex);
  
  double Delta(double x, int d) {
    double pchisq(double q, double df, int lower_tail, int give_log);
    return pchisq(x, (double)d, 1, 0);
  }
 
  double delta(double x, int d) {
    double dchisq(double x, double df, int give_log);
    return dchisq(x, (double)d, 0);
  }
  

  void f1t2(double *z, int n, void *ex2) { // f1v1
    double Delta(double x, int d);
    double delta(double x, int d);
    double x;
    int i, d;    
    x = *((double*)ex2+0);
    d = (int)*((double*)ex2+1);
    for (i = 1; i <= n; i++) {z[i - 1] = Delta(x - z[i - 1], d) * delta(z[i - 1], 1);}
  }


  void f2t1(double *z, int n, void *ex1) { // f1v2
    double Delta(double x, int d);
    double delta(double x, int d);
    int i, nn;    
    nn =  (int)*((double*)ex1+0);
    for (i = 1; i <= n; i++) {z[i - 1] = Delta(2.0 * log(nn) - z[i - 1], 1) * delta(z[i - 1], 1);}
  }
  
  void f2t2(double *z, int n, void *ex2) { // f2v2
    double Delta(double x, int d);
    double delta(double x, int d);
    int i, d;
    double x;
    x = *((double*)ex2+0);
    d = (int)*((double*)ex2+1);
    for (i = 1; i <= n; i++) {z[i - 1] = Delta(x - z[i - 1], 1) * delta(z[i - 1], d);}
  }
  
  void f2t3ab(double *t, int n, void *ex3) { // f3v2
    double Delta(double x, int d);
    double delta(double x, int d);
    int i, d;
    double x, az;
    az = *((double*)ex3+0);
    x = *((double*)ex3+1);
    d = (int)*((double*)ex3+2);
    for (i = 1; i <= n; i++) {t[i - 1] = Delta(x - az - t[i - 1], d) * delta(t[i - 1], 1);}
  }
  
  void f2t3af3t3a(double *z, int n, void *ex3) { // f4v2
    
    void f2t3ab(double *t, int n, void *ex3); // f3v2
    double delta(double x, int d);
    void Rdqags(integr_fn f, void *ex, double *a, double *b,
		double *epsabs, double *epsrel,
		double *result, double *abserr, int *neval, int *ier,
		int *limit, int *lenw, int *last,
		int *iwork, double *work);
    
    int nn, i, d;
    double x;
    double *a, *b, *epsabs, *epsrel, *result, *ex33, *abserr, *work;
    int *last, *limit, *lenw, *ier, *neval, *iwork;
    
    epsabs = new double[1];
    epsrel = new double[1];
    ex33 = new double[3];
    limit = new int[1];
    lenw = new int[1];
    a = new double[1];
    b = new double[1];
    result = new double[1];
    abserr = new double[1];
    last = new int[1];
    ier = new int[1];
    neval = new int[1];

    nn = (int)*((double*)ex3+0);
    x = *((double*)ex3+1);
    d = (int)*((double*)ex3+2);
    
    epsabs[0] = 0.0001220703125;
    epsrel[0] = 0.0001220703125;
    ex33[1] = x;
    ex33[2] = (double)d;
    
    limit[0] = 100;
    lenw[0] = 4 * limit[0];
    
    iwork =  new int[limit[0]];
    work = new double[lenw[0]];
    
    for (i = 1; i <= n; i++) {
      
      a[0] = (2.0 * log(nn) - z[i - 1]);
      b[0] = x - z[i - 1];
      ex33[0] = z[i - 1];
      
      Rdqags(f2t3ab, ex33, a, b, epsabs, epsrel, // f3v2
	     result, abserr, neval, ier,
	     limit, lenw, last,
	     iwork, work);
      
      z[i - 1] = result[0] * delta(z[i - 1], 1);
      
    }
   
    //On libere de la memoire
    delete[] epsabs;
    delete[] epsrel;
    delete[] ex33;
    delete[] limit;
    delete[] lenw;
    delete[] a;
    delete[] b;
    delete[] result;
    delete[] abserr;
    delete[] last;
    delete[] ier;
    delete[] neval;
    delete[] iwork;
    delete[] work;
 
  }
  
  void f2t3bf3t3b(double *z, int n, void *ex3) { // f5v2
    
    void f2t3ab(double *t, int n, void *ex3); // f3v2
    double delta(double x, int d);
    void Rdqags(integr_fn f, void *ex, double *a, double *b,
		double *epsabs, double *epsrel,
		double *result, double *abserr, int *neval, int *ier,
		int *limit, int *lenw, int *last,
		int *iwork, double *work);
    
    int nn, i, d;
    double x;
    double *a, *b, *epsabs, *epsrel, *result, *ex33, *abserr, *work;
    int *last, *limit, *lenw, *ier, *neval, *iwork;
    
    epsabs = new double[1];
    epsrel = new double[1];
    ex33 = new double[3];
    limit = new int[1];
    lenw = new int[1];
    a = new double[1];
    b = new double[1];
    result = new double[1];
    abserr = new double[1];
    last = new int[1];
    ier = new int[1];
    neval = new int[1];

    nn = (int)*((double*)ex3+0);
    x = *((double*)ex3+1);
    d = (int)*((double*)ex3+2);
    
    epsabs[0] = 0.0001220703125;
    epsrel[0] = 0.0001220703125;
    ex33[1] = x;
    ex33[2] = (double)d;
    
    limit[0] = 100;
    lenw[0] = 4 * limit[0];
    a[0] = log(nn);
    
    iwork =   new int[limit[0]];
    work = new double[lenw[0]];

    for (i = 1; i <= n; i++) {
      
      b[0] = x - z[i - 1];
      ex33[0] = z[i - 1];
      
      Rdqags(f2t3ab, ex33, a, b, epsabs, epsrel, // f3v2
	     result, abserr, neval, ier,
	     limit, lenw, last,
	     iwork, work);
      
      z[i-1] = result[0] * delta(z[i - 1], 1);
      
    }
   
    //On libere de la memoire
    delete[] epsabs;
    delete[] epsrel;
    delete[] ex33;
    delete[] limit;
    delete[] lenw;
    delete[] a;
    delete[] b;
    delete[] result;
    delete[] abserr;
    delete[] last;
    delete[] ier;
    delete[] neval;
    delete[] iwork;
    delete[] work;
    
 
  }
  
  void f3t1inside(double *z, int n, void *ex2) { // f1v3
    double Delta(double x, int d);
    double delta(double x, int d);
    double y;
    int i, nn;    
    y = *((double*)ex2+0);
    nn = (int)*((double*)ex2+1);
    for (i = 1; i <= n; i++) {z[i - 1] = Delta(3.0 * log(nn) - y - z[i - 1], 1) * delta(z[i - 1], 1);}
  }


  void f3t1(double *y, int n, void *ex1) { // f2v3
    
    void f3t1inside(double *z, int n, void *ex2); // f1v3
    double delta(double x, int d);
    void Rdqags(integr_fn f, void *ex, double *a, double *b,
		double *epsabs, double *epsrel,
		double *result, double *abserr, int *neval, int *ier,
		int *limit, int *lenw, int *last,
		int *iwork, double *work);
    
    int nn, i;
    double *a, *b, *epsabs, *epsrel, *result, *ex2, *abserr, *work;
    int *last, *limit, *lenw, *ier, *neval, *iwork;
    
    epsabs = new double[1];
    epsrel = new double[1];
    ex2 = new double[3];
    limit = new int[1];
    lenw = new int[1];
    a = new double[1];
    b = new double[1];
    result = new double[1];
    abserr = new double[1];
    last = new int[1];
    ier = new int[1];
    neval = new int[1];

    nn = (int)*((double*)ex1+0);
    
    epsabs[0] = 0.0001220703125;
    epsrel[0] = 0.0001220703125;
    
    limit[0] = 100;
    lenw[0] = 4 * limit[0];
    
    iwork =  new int[*(limit+0)];
    work = new double[*(lenw+0)];
    
    ex2[1] = (double)nn;
    
    for (i = 1; i <= n; i++) {
      
      a[0] = 0;
      b[0] = 2.0 * log(nn) - y[i - 1];
      ex2[0] = y[i - 1];
      
      Rdqags(f3t1inside, ex2, a, b, epsabs, epsrel, // f1v3
	     result, abserr, neval, ier,
	     limit, lenw, last,
	     iwork, work);
      
      y[i - 1] = result[0] * delta(y[i - 1], 1);
      
    }
   
    //On libere de la memoire
    delete[] epsabs;
    delete[] epsrel;
    delete[] ex2;
    delete[] limit;
    delete[] lenw;
    delete[] a;
    delete[] b;
    delete[] result;
    delete[] abserr;
    delete[] last;
    delete[] ier;
    delete[] neval;
    delete[] iwork;
    delete[] work;
 
  }


  void f3t2a(double *z, int n, void *ex3) { // f3v3
    double Delta(double x, int d);
    double delta(double x, int d);
    double x;
    int i, d; // nn;    
    //    nn = (int)*((double*)ex3+0);
    x = *((double*)ex3+1);
    d = (int)*((double*)ex3+2);
    for (i = 1; i <= n; i++) {z[i - 1] = Delta(x - z[i - 1], d) * delta(z[i - 1], 1);}
  }

  void f3t2b(double *z, int n, void *ex1) { // f4v3
    double Delta(double x, int d);
    double delta(double x, int d);
    int i, nn;    
    nn = (int)*((double*)ex1+0);
    for (i = 1; i <= n; i++) {z[i - 1] = Delta(2.0 * log(nn) - z[i - 1], 1) * delta(z[i - 1], 1);}
  }

  void f3t4a(double *x1, int n, void *ex3) {
    
    void f3t4ainside1(double *x2, int n, void *ex4);
    double delta(double x, int d);
    void Rdqags(integr_fn f, void *ex, double *a, double *b,
		double *epsabs, double *epsrel,
		double *result, double *abserr, int *neval, int *ier,
		int *limit, int *lenw, int *last,
		int *iwork, double *work);
    
    int nn, i, d;
    double x;
    double *a, *b, *epsabs, *epsrel, *result, *ex41, *abserr, *work;
    int *last, *limit, *lenw, *ier, *neval, *iwork;
    
    epsabs = new double[1];
    epsrel = new double[1];
    ex41 = new double[4];
    limit = new int[1];
    lenw = new int[1];
    a = new double[1];
    b = new double[1];
    result = new double[1];
    abserr = new double[1];
    last = new int[1];
    ier = new int[1];
    neval = new int[1];

    nn = (int)*((double*)ex3+0);
    x = *((double*)ex3+1);
    d = (int)*((double*)ex3+2);
    
    epsabs[0] = 0.0001220703125;
    epsrel[0] = 0.0001220703125;
    ex41[0] = (double)nn;
    ex41[1] = x;
    ex41[2] = (double)d;
    
    limit[0] = 100;
    lenw[0] = 4 * limit[0];
    
    iwork =   new int[limit[0]];
    work = new double[lenw[0]];

    a[0] = 0.0;
    b[0] = log(nn);

    for (i = 1; i <= n; i++) {
      
      ex41[3] = x1[i - 1];
      
      Rdqags(f3t4ainside1, ex41, a, b, epsabs, epsrel,
	     result, abserr, neval, ier,
	     limit, lenw, last,
	     iwork, work);
      
      x1[i-1] = result[0] * delta(x1[i - 1], 1);
      
    }
   
    //On libere de la memoire
    delete[] epsabs;
    delete[] epsrel;
    delete[] ex41;
    delete[] limit;
    delete[] lenw;
    delete[] a;
    delete[] b;
    delete[] result;
    delete[] abserr;
    delete[] last;
    delete[] ier;
    delete[] neval;
    delete[] iwork;
    delete[] work;
    
 
  }

  void f3t4ainside1(double *x2, int n, void *ex4) {
    
    void f3t4ainside2(double *x3, int n, void *ex4);
    double delta(double x, int d);
    void Rdqags(integr_fn f, void *ex, double *a, double *b,
		double *epsabs, double *epsrel,
		double *result, double *abserr, int *neval, int *ier,
		int *limit, int *lenw, int *last,
		int *iwork, double *work);
    
    int nn, i, d;
    double x, x1;
    double *a, *b, *epsabs, *epsrel, *result, *ex42, *abserr, *work;
    int *last, *limit, *lenw, *ier, *neval, *iwork;
    
    epsabs = new double[1];
    epsrel = new double[1];
    ex42 = new double[4];
    limit = new int[1];
    lenw = new int[1];
    a = new double[1];
    b = new double[1];
    result = new double[1];
    abserr = new double[1];
    last = new int[1];
    ier = new int[1];
    neval = new int[1];

    nn = (int)*((double*)ex4+0);
    x = *((double*)ex4+1);
    d = (int)*((double*)ex4+2);
    x1 = *((double*)ex4+3);
    
    epsabs[0] = 0.0001220703125;
    epsrel[0] = 0.0001220703125;
    ex42[1] = x;
    ex42[2] = (double)d;
    ex42[3] = x1;
    
    limit[0] = 100;
    lenw[0] = 4 * limit[0];
    
    iwork =   new int[limit[0]];
    work = new double[lenw[0]];

    for (i = 1; i <= n; i++) {
      
      a[0] = 3.0 * log(nn) - x1 - x2[i - 1];
      b[0] = x - x1 - x2[i - 1];
      ex42[0] = x2[i - 1];
      
      Rdqags(f3t4ainside2, ex42, a, b, epsabs, epsrel,
	     result, abserr, neval, ier,
	     limit, lenw, last,
	     iwork, work);
      
      x2[i-1] = result[0] * delta(x2[i - 1], 1);
      
    }
   
    //On libere de la memoire
    delete[] epsabs;
    delete[] epsrel;
    delete[] ex42;
    delete[] limit;
    delete[] lenw;
    delete[] a;
    delete[] b;
    delete[] result;
    delete[] abserr;
    delete[] last;
    delete[] ier;
    delete[] neval;
    delete[] iwork;
    delete[] work;
    
 
  }

  void f3t4ainside2(double *x3, int n, void *ex4) {
    double Delta(double x, int d);
    double delta(double x, int d);
    int i, d;
    double x, x1, x2;    
    x2 = *((double*)ex4+0);
    x = *((double*)ex4+1);
    d = (int)*((double*)ex4+2);
    x1 = *((double*)ex4+3);
    for (i = 1; i <= n; i++) {x3[i - 1] = Delta(x - x1 - x2 - x3[i - 1], d) * delta(x3[i - 1], 1);}
  }

  void f3t4b(double *x1, int n, void *ex3) {
    
    void f3t4binside1(double *x2, int n, void *ex4);
    double delta(double x, int d);
    void Rdqags(integr_fn f, void *ex, double *a, double *b,
		double *epsabs, double *epsrel,
		double *result, double *abserr, int *neval, int *ier,
		int *limit, int *lenw, int *last,
		int *iwork, double *work);
    
    int nn, i, d;
    double x;
    double *a, *b, *epsabs, *epsrel, *result, *ex43, *abserr, *work;
    int *last, *limit, *lenw, *ier, *neval, *iwork;
    
    epsabs = new double[1];
    epsrel = new double[1];
    ex43 = new double[4];
    limit = new int[1];
    lenw = new int[1];
    a = new double[1];
    b = new double[1];
    result = new double[1];
    abserr = new double[1];
    last = new int[1];
    ier = new int[1];
    neval = new int[1];

    nn = (int)*((double*)ex3+0);
    x = *((double*)ex3+1);
    d = (int)*((double*)ex3+2);
    
    epsabs[0] = 0.0001220703125;
    epsrel[0] = 0.0001220703125;
    ex43[0] = (double)nn;
    ex43[1] = x;
    ex43[2] = (double)d;
    
    limit[0] = 100;
    lenw[0] = 4 * limit[0];
    
    iwork =   new int[limit[0]];
    work = new double[lenw[0]];

    for (i = 1; i <= n; i++) {
      
      a[0] = 2.0 * log(nn) - x1[i - 1];
      b[0] = x - log(nn) - x1[i - 1];
      ex43[3] = x1[i - 1];
      
      Rdqags(f3t4binside1, ex43, a, b, epsabs, epsrel,
	     result, abserr, neval, ier,
	     limit, lenw, last,
	     iwork, work);
      
      x1[i-1] = result[0] * delta(x1[i - 1], 1);
      
    }
   
    //On libere de la memoire
    delete[] epsabs;
    delete[] epsrel;
    delete[] ex43;
    delete[] limit;
    delete[] lenw;
    delete[] a;
    delete[] b;
    delete[] result;
    delete[] abserr;
    delete[] last;
    delete[] ier;
    delete[] neval;
    delete[] iwork;
    delete[] work;
    
 
  }

  void f3t4binside1(double *x2, int n, void *ex4) {
    
    void f3t4binside2(double *x3, int n, void *ex4);
    double delta(double x, int d);
    void Rdqags(integr_fn f, void *ex, double *a, double *b,
		double *epsabs, double *epsrel,
		double *result, double *abserr, int *neval, int *ier,
		int *limit, int *lenw, int *last,
		int *iwork, double *work);
    
    int nn, i, d;
    double x, x1;
    double *a, *b, *epsabs, *epsrel, *result, *ex44, *abserr, *work;
    int *last, *limit, *lenw, *ier, *neval, *iwork;
    
    epsabs = new double[1];
    epsrel = new double[1];
    ex44 = new double[4];
    limit = new int[1];
    lenw = new int[1];
    a = new double[1];
    b = new double[1];
    result = new double[1];
    abserr = new double[1];
    last = new int[1];
    ier = new int[1];
    neval = new int[1];

    nn = (int)*((double*)ex4+0);
    x = *((double*)ex4+1);
    d = (int)*((double*)ex4+2);
    x1 = *((double*)ex4+3);
    
    epsabs[0] = 0.0001220703125;
    epsrel[0] = 0.0001220703125;
    ex44[1] = x;
    ex44[2] = (double)d;
    ex44[3] = x1;
    
    limit[0] = 100;
    lenw[0] = 4 * limit[0];
    
    iwork =   new int[limit[0]];
    work = new double[lenw[0]];

    for (i = 1; i <= n; i++) {
      
      a[0] = log(nn);
      b[0] = x - x1 - x2[i - 1];
      ex44[0] = x2[i - 1];
      
      Rdqags(f3t4ainside2, ex44, a, b, epsabs, epsrel,
	     result, abserr, neval, ier,
	     limit, lenw, last,
	     iwork, work);
      
      x2[i-1] = result[0] * delta(x2[i - 1], 1);
      
    }
   
    //On libere de la memoire
    delete[] epsabs;
    delete[] epsrel;
    delete[] ex44;
    delete[] limit;
    delete[] lenw;
    delete[] a;
    delete[] b;
    delete[] result;
    delete[] abserr;
    delete[] last;
    delete[] ier;
    delete[] neval;
    delete[] iwork;
    delete[] work;
    
 
  }

  void f3t4c(double *x1, int n, void *ex3) {
    
    void f3t4ainside1(double *x2, int n, void *ex4);
    double delta(double x, int d);
    void Rdqags(integr_fn f, void *ex, double *a, double *b,
		double *epsabs, double *epsrel,
		double *result, double *abserr, int *neval, int *ier,
		int *limit, int *lenw, int *last,
		int *iwork, double *work);
    
    int nn, i, d;
    double x;
    double *a, *b, *epsabs, *epsrel, *result, *ex41, *abserr, *work;
    int *last, *limit, *lenw, *ier, *neval, *iwork;
    
    epsabs = new double[1];
    epsrel = new double[1];
    ex41 = new double[4];
    limit = new int[1];
    lenw = new int[1];
    a = new double[1];
    b = new double[1];
    result = new double[1];
    abserr = new double[1];
    last = new int[1];
    ier = new int[1];
    neval = new int[1];

    nn = (int)*((double*)ex3+0);
    x = *((double*)ex3+1);
    d = (int)*((double*)ex3+2);
    
    epsabs[0] = 0.0001220703125;
    epsrel[0] = 0.0001220703125;
    ex41[0] = (double)nn;
    ex41[1] = x;
    ex41[2] = (double)d;
    
    limit[0] = 100;
    lenw[0] = 4 * limit[0];
    
    iwork =   new int[limit[0]];
    work = new double[lenw[0]];

    a[0] = log(nn);

    for (i = 1; i <= n; i++) {
      
      b[0] = 2.0 * log(nn) - x1[i - 1];

      ex41[3] = x1[i - 1];
      
      Rdqags(f3t4ainside1, ex41, a, b, epsabs, epsrel,
	     result, abserr, neval, ier,
	     limit, lenw, last,
	     iwork, work);
      
      x1[i-1] = result[0] * delta(x1[i - 1], 1);
      
    }
   
    //On libere de la memoire
    delete[] epsabs;
    delete[] epsrel;
    delete[] ex41;
    delete[] limit;
    delete[] lenw;
    delete[] a;
    delete[] b;
    delete[] result;
    delete[] abserr;
    delete[] last;
    delete[] ier;
    delete[] neval;
    delete[] iwork;
    delete[] work;
    
 
  }

  void f3t4d(double *x1, int n, void *ex3) {
    
    void f3t4cinside1(double *x2, int n, void *ex4);
    double delta(double x, int d);
    void Rdqags(integr_fn f, void *ex, double *a, double *b,
		double *epsabs, double *epsrel,
		double *result, double *abserr, int *neval, int *ier,
		int *limit, int *lenw, int *last,
		int *iwork, double *work);
    
    int nn, i, d;
    double x;
    double *a, *b, *epsabs, *epsrel, *result, *ex41, *abserr, *work;
    int *last, *limit, *lenw, *ier, *neval, *iwork;
    
    epsabs = new double[1];
    epsrel = new double[1];
    ex41 = new double[4];
    limit = new int[1];
    lenw = new int[1];
    a = new double[1];
    b = new double[1];
    result = new double[1];
    abserr = new double[1];
    last = new int[1];
    ier = new int[1];
    neval = new int[1];

    nn = (int)*((double*)ex3+0);
    x = *((double*)ex3+1);
    d = (int)*((double*)ex3+2);
    
    epsabs[0] = 0.0001220703125;
    epsrel[0] = 0.0001220703125;
    ex41[0] = (double)nn;
    ex41[1] = x;
    ex41[2] = (double)d;
    
    limit[0] = 100;
    lenw[0] = 4 * limit[0];
    
    iwork =   new int[limit[0]];
    work = new double[lenw[0]];

    a[0] = 0.0;
    b[0] = log(nn);

    for (i = 1; i <= n; i++) {
      
      ex41[3] = x1[i - 1];
      
      Rdqags(f3t4cinside1, ex41, a, b, epsabs, epsrel,
	     result, abserr, neval, ier,
	     limit, lenw, last,
	     iwork, work);
      
      x1[i-1] = result[0] * delta(x1[i - 1], 1);
      
    }
   
    //On libere de la memoire
    delete[] epsabs;
    delete[] epsrel;
    delete[] ex41;
    delete[] limit;
    delete[] lenw;
    delete[] a;
    delete[] b;
    delete[] result;
    delete[] abserr;
    delete[] last;
    delete[] ier;
    delete[] neval;
    delete[] iwork;
    delete[] work;
    
 
  }

  void f3t4cinside1(double *x2, int n, void *ex4) {
    
    void f3t4ainside2(double *x3, int n, void *ex4);
    double delta(double x, int d);
    void Rdqags(integr_fn f, void *ex, double *a, double *b,
		double *epsabs, double *epsrel,
		double *result, double *abserr, int *neval, int *ier,
		int *limit, int *lenw, int *last,
		int *iwork, double *work);
    
    int nn, i, d;
    double x, x1;
    double *a, *b, *epsabs, *epsrel, *result, *ex42, *abserr, *work;
    int *last, *limit, *lenw, *ier, *neval, *iwork;
    
    epsabs = new double[1];
    epsrel = new double[1];
    ex42 = new double[4];
    limit = new int[1];
    lenw = new int[1];
    a = new double[1];
    b = new double[1];
    result = new double[1];
    abserr = new double[1];
    last = new int[1];
    ier = new int[1];
    neval = new int[1];

    nn = (int)*((double*)ex4+0);
    x = *((double*)ex4+1);
    d = (int)*((double*)ex4+2);
    x1 = *((double*)ex4+3);
    
    epsabs[0] = 0.0001220703125;
    epsrel[0] = 0.0001220703125;
    ex42[1] = x;
    ex42[2] = (double)d;
    ex42[3] = x1;
    
    limit[0] = 100;
    lenw[0] = 4 * limit[0];
    
    iwork =   new int[limit[0]];
    work = new double[lenw[0]];

    for (i = 1; i <= n; i++) {
      
      a[0] = 2.0 * log(nn) - x2[i - 1];
      b[0] = x - x1 - x2[i - 1];
      ex42[0] = x2[i - 1];
      
      Rdqags(f3t4ainside2, ex42, a, b, epsabs, epsrel,
	     result, abserr, neval, ier,
	     limit, lenw, last,
	     iwork, work);
      
      x2[i-1] = result[0] * delta(x2[i - 1], 1);
      
    }
   
    //On libere de la memoire
    delete[] epsabs;
    delete[] epsrel;
    delete[] ex42;
    delete[] limit;
    delete[] lenw;
    delete[] a;
    delete[] b;
    delete[] result;
    delete[] abserr;
    delete[] last;
    delete[] ier;
    delete[] neval;
    delete[] iwork;
    delete[] work;
    
 
  }

  void f3t4e(double *x1, int n, void *ex3) {
    
    void f3t4binside1(double *x2, int n, void *ex4);
    double delta(double x, int d);
    void Rdqags(integr_fn f, void *ex, double *a, double *b,
		double *epsabs, double *epsrel,
		double *result, double *abserr, int *neval, int *ier,
		int *limit, int *lenw, int *last,
		int *iwork, double *work);
    
    int nn, i, d;
    double x;
    double *a, *b, *epsabs, *epsrel, *result, *ex43, *abserr, *work;
    int *last, *limit, *lenw, *ier, *neval, *iwork;
    
    epsabs = new double[1];
    epsrel = new double[1];
    ex43 = new double[4];
    limit = new int[1];
    lenw = new int[1];
    a = new double[1];
    b = new double[1];
    result = new double[1];
    abserr = new double[1];
    last = new int[1];
    ier = new int[1];
    neval = new int[1];

    nn = (int)*((double*)ex3+0);
    x = *((double*)ex3+1);
    d = (int)*((double*)ex3+2);
    
    epsabs[0] = 0.0001220703125;
    epsrel[0] = 0.0001220703125;
    ex43[0] = (double)nn;
    ex43[1] = x;
    ex43[2] = (double)d;
    
    limit[0] = 100;
    lenw[0] = 4 * limit[0];
    
    iwork =   new int[limit[0]];
    work = new double[lenw[0]];

    a[0] = log(nn);

    for (i = 1; i <= n; i++) {
      
      b[0] = x - log(nn) - x1[i - 1];
      ex43[3] = x1[i - 1];
      
      Rdqags(f3t4binside1, ex43, a, b, epsabs, epsrel,
	     result, abserr, neval, ier,
	     limit, lenw, last,
	     iwork, work);
      
      x1[i-1] = result[0] * delta(x1[i - 1], 1);
      
    }
   
    //On libere de la memoire
    delete[] epsabs;
    delete[] epsrel;
    delete[] ex43;
    delete[] limit;
    delete[] lenw;
    delete[] a;
    delete[] b;
    delete[] result;
    delete[] abserr;
    delete[] last;
    delete[] ier;
    delete[] neval;
    delete[] iwork;
    delete[] work;
    
 
  }





  void R_cumsum(double *x, int *n, double *na_value, double *ans) {
    double sum;
    int i;

    for(i = 0; i < *n; i++) ans[i] = *na_value;

    sum = 0.0;
    for(i = 0; i < *n; i++) {
      if(x[i] == *na_value) break;
      sum += x[i];
      ans[i] = sum;
    }
  } 
  
} // extern C

