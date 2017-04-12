//#include<math.h>

void H1(int n, double *x) {
  int i;
  for (i=0;i<n;i++) {
    x[i] = sqrt(3.0) * x[i];
  }
  return;
}

