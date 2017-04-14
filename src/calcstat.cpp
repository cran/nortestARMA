// g++ -c calcstat.cpp -o calcstat.o -I"/usr/lib/R/include"
// g++ -shared -o calcstatC.so calcstat.o -I"/usr/lib/R/include" -L"/usr/lib" -lR

#include <R.h>
#include "Rmath.h"
//#include <R_ext/Applic.h>

using std::sqrt;
using std::log;
using std::pow;

extern "C" {
    

#include "polys/H1.cpp"
#include "polys/H2.cpp"
#include "polys/H3.cpp"
#include "polys/H4.cpp"
#include "polys/H5.cpp"
#include "polys/H6.cpp"
#include "polys/H7.cpp"
#include "polys/H8.cpp"
#include "polys/H9.cpp"
#include "polys/H10.cpp"
  
  /*
#include "H1etoile.cc"
#include "H2etoile.cc"
#include "H3etoile.cc"
#include "H4etoile.cc"
#include "H5etoile.cc"
#include "H6etoile.cc"
#include "H7etoile.cc"
#include "H8etoile.cc"
#include "H9etoile.cc"
#include "H10etoile.cc"
  */

  //#include "quantile.cpp"

  // #include <stdio.h>
  
  //  void F77_SUB(calcstat)(double *residu, double *sigmahat, int *nT, int *d, double *RK, double *RKduch, double *RKtagne, int *Kopt, int *Koptduch, int *Kopttagne, double *st, double *stduch, double *sttagne, double *pval, double *pvalduch, double *pvaltagne)

  void calcstat(double *residu, double *sigmahat, int *nT, int *d, double *RK, double *RKduch, double *RKtagne, int *Kopt, int *Koptduch, int *Kopttagne, double *st, double *stduch, double *sttagne, double *pval, double *pvalduch, double *pvaltagne)

{

  double pnorm(double x, double mu, double sigma, int lower_tail,int give_log);
  double pvalst(double* x, int *xlen, int* n, int* d, int* K, double *res, double* term1, double* term2, double* term3, double* term4);
  void R_cumsum(double *x, int *n, double *na_value, double *ans);

  void H1(int n, double *U);
  void H2(int n, double *U);
  void H3(int n, double *U);
  void H4(int n, double *U);
  void H5(int n, double *U);
  void H6(int n, double *U);
  void H7(int n, double *U);
  void H8(int n, double *U);
  void H9(int n, double *U);
  void H10(int n, double *U);
  /*
  void H1etoi(int n, double *U);
  void H2etoi(int n, double *U);
  void H3etoi(int n, double *U);
  void H4etoi(int n, double *U);
  void H5etoi(int n, double *U);
  void H6etoi(int n, double *U);
  void H7etoi(int n, double *U);
  void H8etoi(int n, double *U);
  void H9etoi(int n, double *U);
  void H10eto(int n, double *U);
  */

  int *BigD, k, t, n, K, *xlen, *Kpvalst;
  xlen = new int[1];
  xlen[0] = (int)1;
  Kpvalst = new int[1];

  double *res, *term1, *term2, *term3, *term4;
  res = new double[1];
  term1 = new double[1];
  term2 = new double[1];
  term3 = new double[1];
  term4 = new double[1];
  res[0] = 0.0;
  term1[0] = 0.0; 
  term2[0] = 0.0; 
  term3[0] = 0.0; 
  term4[0] = 0.0; 

  n = nT[0];

  //  BigD=(int)log(double(n));
  BigD = new int[1];
  BigD[0] = 10;

  double rktemp1, rktemp2, rktemp3, penal, mxLed, mxLedduch, mxLedtagne; // ,rktemp, maxbn; //, c;

  //  c = 2.4;
  //  penal = (1.0+(std::sqrt(2.0)*exp(0)/double(n)))*log(double(n));
  // penal = log(double(n));
  // penal = 2;
	
  double *U, *Uavant, *hU, *Ledwi, *Ledwiduch, *Ledwitagne, *hKtemp; // *hUetoi, 
  U = new double[n];
  Uavant = new double[n];
  //  hUetoi = new double[n * BigD[0]];
  hU = new double[n*BigD[0]];
  Ledwi = new double[BigD[0]];
  Ledwiduch = new double[BigD[0]];
  Ledwitagne = new double[BigD[0]];
  hKtemp = new double[BigD[0]];

  double *a, *b, *a2, *b2, *cumsuma, *cumsuma2, *cumsumb, *cumsumb2, *penalab, *na_value;
  a = new double[BigD[0]];
  a2 = new double[BigD[0]];
  b = new double[BigD[0]];
  b2 = new double[BigD[0]];
  cumsuma = new double[BigD[0]];
  cumsuma2 = new double[BigD[0]];
  cumsumb = new double[BigD[0]];
  cumsumb2 = new double[BigD[0]];
  penalab = new double[BigD[0]];
  na_value = new double[BigD[0]];

  a[0] = 0.9772050238058398; // C'est sqrt(3.0/pi). Avant j'avais mis: 0.9772049160826732;
  a[1] = 0.0;
  a[2] = 0.18300807767356908;
  a[3] = 0.0;
  a[4] = 0.08169877690412206;
  a[5] = 0.0;
  a[6] = 0.04772914572835635;
  a[7] = 0.0;
  a[8] = 0.03188018904558847;
  a[9] = 0.0;

  b[0] = 0.0;
  b[1] = 1.232808888123; // C'est sqrt(15.0)/pi. Avant, j'avais mis: 1.2328087469140054;
  b[2] = 0.0;
  b[3] = 0.5211244103155082;
  b[4] = 0.0;
  b[5] = 0.3045142864019812;
  b[6] = 0.0;
  b[7] = 0.20558879812557124;
  b[8] = 0.0;
  b[9] = 0.1507705607248925;

  for (k = 1; k <= BigD[0]; k++) na_value[k-1] = -1.0;

  for (k = 1; k <= BigD[0]; k++) a2[k-1] = pow(a[k-1], 2.0);
  for (k = 1; k <= BigD[0]; k++) b2[k-1] = pow(b[k-1], 2.0);
	    
  R_cumsum(a, BigD, na_value, cumsuma);
  R_cumsum(a2, BigD, na_value, cumsuma2);
  R_cumsum(b, BigD, na_value, cumsumb);
  R_cumsum(b2, BigD, na_value, cumsumb2);

  for (k = 1; k <= BigD[0]; k++) penalab[k-1] = cumsumb[k-1] / std::sqrt(2.0 - cumsumb2[k-1]) + cumsuma[k-1] / std::sqrt(1.0 - cumsuma2[k-1]);


  for (t = 1; t <= n; t++) {U[t-1] = 2.0 *  pnorm(residu[t-1] / sigmahat[0], 0.0, 1.0, 1,0) - 1.0;}

  /*

  for (k = 1; k <= BigD[0]; k++) {
  
    for (t = 1; t <= n; t++) {Uavant[t - 1] = U[t - 1]; }

    if (k == 1)  H1etoi(n, U);
    if (k == 2)  H2etoi(n, U);
    if (k == 3)  H3etoi(n, U);
    if (k == 4)  H4etoi(n, U);
    if (k == 5)  H5etoi(n, U);
    if (k == 6)  H6etoi(n, U);
    if (k == 7)  H7etoi(n, U);
    if (k == 8)  H8etoi(n, U);
    if (k == 9)  H9etoi(n, U);
    if (k == 10) H10eto(n, U);
		
    for (t = 1; t <= n; t++) {hUetoi[t - 1 + (k - 1) * n] = U[t - 1];}

    for (t = 1; t <= n; t++) {U[t - 1] = Uavant[t - 1]; }

  }

  */


  for (k = 1; k <= BigD[0]; k++) {
  
    for (t = 1; t <= n; t++) Uavant[t - 1] = U[t - 1]; 

    if (k == 1)  H1(n, U);
    if (k == 2)  H2(n, U);
    if (k == 3)  H3(n, U);
    if (k == 4)  H4(n, U);
    if (k == 5)  H5(n, U);
    if (k == 6)  H6(n, U);
    if (k == 7)  H7(n, U);
    if (k == 8)  H8(n, U);
    if (k == 9)  H9(n, U);
    if (k == 10) H10(n, U);
		
    for (t = 1; t <= n; t++) hU[t - 1 + (k - 1) * n] = U[t - 1];

    for (t = 1; t <= n; t++) U[t - 1] = Uavant[t - 1];

  }
  
    /*

  //      RKtagne[k]:c'est la statistique de mon test calculee avec les hketoile
  rktemp = 0.0;

  for (k = 1; k <= BigD[0]; k++) {

    hKtemp[k - 1] = 0.0;

    for (t = 1; t <= n; t++) {hKtemp[k - 1] = hKtemp[k - 1] + hUetoi[t - 1 + (k - 1) * n];}

    rktemp = rktemp + (pow(hKtemp[k - 1], 2.0)) / (double(n));
    RKtagne[k - 1] = rktemp;
        printf("%8.4f \n", rktemp);

  }


  //      RK[k]:c'est la statistique de mon test calculee avec les hk
  rktemp = 0.0;

  for (k = 1; k <= BigD[0]; k++) {

    hKtemp[k-1] = 0.0;

    for (t = 1; t <= n; t++) {hKtemp[k - 1] = hKtemp[k - 1] + hU[t - 1 + (k - 1) * n];}

    rktemp = rktemp + (pow(hKtemp[k - 1], 2.0)) / (double(n));
    RK[k - 1] = rktemp;

    //    printf("%8.2f \n", rktemp);

  }
  

  //      On essaye en utilisant l'equation (22) du nouveau papier avec Duchesne

  for (k = 1; k <= BigD[0]; k++) {

    hKtemp[k - 1] = 0.0;

    for (t = 1; t <= n; t++) {hKtemp[k - 1] = hKtemp[k - 1] + hU[t - 1 + (k - 1) * n];}
    //      printf("%8.8f \n", hKtemp[k - 1]);

  }
  */



  for (K = 1; K <= BigD[0]; K++) {
    hKtemp[K-1] = 0.0;
    for (t = 1; t <= n; t++) {hKtemp[K - 1] = hKtemp[K - 1] + hU[t - 1 + (K - 1) * n];}

    rktemp1 = 0.0;
    rktemp2 = 0.0;
    rktemp3 = 0.0;
    for (k = 1; k <= K; k++) {
      rktemp1 = rktemp1 + pow(hKtemp[k - 1], 2.0);
      rktemp2 = rktemp2 + a[k - 1] * hKtemp[k - 1];
      rktemp3 = rktemp3 + b[k - 1] * hKtemp[k - 1];
    }
    RK[K - 1] = rktemp1 / (double(n)) ;
    RKduch[K - 1] = (rktemp1 + pow(rktemp3, 2.0) / (2.0 - cumsumb2[K - 1])) / (double(n)) ;
    RKtagne[K - 1] = (rktemp1 + pow(rktemp2, 2.0) / (1.0 - cumsuma2[K - 1]) + pow(rktemp3, 2.0) / (2.0 - cumsumb2[K - 1])) / (double(n)) ;

    //      printf("%8.4f \n", Rtagne[K - 1]);

  }



  //  hKtemp[j-1]/std::sqrt(double(n)) : std::sqrt(n)*bhatnj comme dans Inglot (2006)
  //   maxbn = fabs(hKtemp[0]/std::sqrt(double(n)));

  // for (k=*(d+0);k<=BigD[0]-1;k++) {
  //   if (fabs(hKtemp[k]/std::sqrt(double(n))) > maxbn) maxbn = fabs(hKtemp[k]/std::sqrt(double(n)));
  //  }

  //  if (maxbn <= std::sqrt(c*log(double(n)))) {penal = log(double(n));} else penal = 2.0;


  penal = std::log(double(n));

  //      On calcule le K optimal par la methode de Ledwina: (critere S2)
     //   for (k=*(d+0);k<=BigD[0];k++) {Ledwitagne[k-1]=RKtagne[k-1]-(k+penalab[k-1])*penal;}
  for (k = d[0]; k <= BigD[0]; k++) {
    Ledwi[k - 1] = RK[k - 1] - k * penal;
    Ledwiduch[k - 1] = RKduch[k - 1] - k * penal;
    Ledwitagne[k - 1] = RKtagne[k - 1] - k * penal;
  }
  mxLed = Ledwi[d[0] - 1];
  mxLedduch = Ledwiduch[d[0] - 1];
  mxLedtagne = Ledwitagne[d[0] - 1];
  Kopt[0] = d[0];
  Koptduch[0] = d[0];
  Kopttagne[0] = d[0];
  for (k = d[0]; k <= BigD[0] - 1; k++) {
    if (Ledwi[k] > mxLed) {
      mxLed = Ledwi[k];
      Kopt[0] = k + 1;
    }
    if (Ledwiduch[k] > mxLedduch) {
      mxLedduch = Ledwiduch[k];
      Koptduch[0] = k + 1;
    }
    if (Ledwitagne[k] > mxLedtagne) {
      mxLedtagne = Ledwitagne[k];
      Kopttagne[0] = k + 1;
    }
  }
  
  st[0] = RK[Kopt[0] - 1];
  stduch[0] = RKduch[Koptduch[0] - 1];
  sttagne[0] = RKtagne[Kopttagne[0] - 1];


  // pvaleurstat
  Kpvalst[0] = d[0] - 1;
  pvalst(st, xlen, nT, d, Kpvalst, res, term1, term2, term3, term4);
  pval[0] = res[0];

  Kpvalst[0] = d[0];
  pvalst(stduch, xlen, nT, d, Kpvalst, res, term1, term2, term3, term4);
  pvalduch[0] = res[0];

  Kpvalst[0] = Kopttagne[0];
  pvalst(sttagne, xlen, nT, d, Kpvalst, res, term1, term2, term3, term4);
  pvaltagne[0] = res[0];



  //On libere de la memoire
  delete[] U;
  delete[] Uavant;
  //  delete[] hUetoi;
  delete[] Ledwi;
  delete[] Ledwiduch;
  delete[] Ledwitagne;
  delete[] hU;
  delete[] hKtemp;

  delete[] BigD;
  delete[] a;
  delete[] a2;
  delete[] b;
  delete[] b2;
  delete[] cumsuma;
  delete[] cumsuma2;
  delete[] cumsumb;
  delete[] cumsumb2;
  delete[] penalab;
  delete[] na_value;

  delete[] xlen;
  delete[] Kpvalst;
  delete[] res;
  delete[] term1;
  delete[] term2;
  delete[] term3;
  delete[] term4;


}

  
  
  
} // extern C

