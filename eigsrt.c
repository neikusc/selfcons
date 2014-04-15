#include "math.h"
#define NRANSI
// #include "nrutil.h"

/* sort eigen values in decending order */
void eigsrt(double d[], double **v, double n) {
  int k,j,i;
  double p;
  
  for (i=1;i<n;i++) {
    p=d[k=i];
    for (j=i+1;j<=n;j++)
      if (d[j] >= p) p=d[k=j];
    if (k != i) {
      d[k]=d[i];
      d[i]=p;
      for (j=1;j<=n;j++) {
	p=v[j][i];
	v[j][i]=v[j][k];
	v[j][k]=p;
      }
    }
  }
  
}

/*-------------------------------------------------------
  Include a matrix and its inverse
-------------------------------------------------------*/

void eigvec(double **h, double **v, double **vt, int n) {
  int i, j;

  for(i=1; i<=n; i++)
    for(j=1; j<=n; j++) {
      v[i][j] = h[i][j];
      vt[j][i] = h[i][j];
    }
}

/*-------------------------------------------------------
  Product of two matrices
-------------------------------------------------------*/
/*
void MatMulti(double **v, double **vt, double **m, int n) {
  int i, j, k;
  double sum;
  for(i=1; i<=n; i++)
  for(j=1; j<=n; j++){
    sum =0;
    for(k=1; k<=n; k++) {
      sum = sum + vt[i][k]*v[k][j];
    }
    m[i][j] = sum;
  }
}
*/
