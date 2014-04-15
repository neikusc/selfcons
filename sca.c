
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "nrutil.c"
#include "nrutil.h"
#include "tred2.c"
#include "tqli.c"
#include "eigsrt.c"
#include "shell.c"
#include "sca.h"
#include "moreo.in"

int main(int argc, char **argv) {

  FILE *p1;
  p1 = fopen("para.dat","w");

  h  = dmatrix(1,N,1,N);
  V  = dmatrix(1,N,1,N);
  Vt = dmatrix(1,N,1,N);
  N1 = dmatrix(1,N,1,N);
  N2 = dmatrix(1,N,1,N);
  M1 = dmatrix(1,N,1,N);
  M2 = dmatrix(1,N,1,N);
  d  = dvector(1,N);
  e  = dvector(1,N);

  MPI_Init(&argc,&argv); /* Initialize the MPI environment */
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);  /* My processor ID */
  MPI_Comm_size(MPI_COMM_WORLD, &nproc); /* Return the number of proccesors */

  /* Vector index of processors */
  vid[0] = myid%vproc[1];
  vid[1] = myid/vproc[1];
  Lx = pi/vproc[0];   
  Ly = 2*pi/vproc[1]; 
  kx = Lx/Nx;
  ky = Ly/Ny;

  U = 2;
  J = U/4;
  Uh = (U-2.5*J);
  
  CompPara();

  printf("\n Hey Kien, It is done ...");
  return 0;
}
/*-----------------------------------------------------------
  Computing parameters n1, n2, m1, m2
-----------------------------------------------------------*/
void CompPara() {
  n = 0.5;
  n1old = 0.5;
  n2old = 0.5;
  m1old = 0.8;
  m2old = 0.8;

  Nmax = 30;
  for(sweep = 1; sweep <= Nmax; sweep++) {
    n1new = 0;
    n2new = 0;
    m1new = 0;
    m2new = 0;
    eSearch();

    for(ix = 0; ix < Nx; ix++)
    for(iy = 0; iy < Ny; iy++){
      px = -pi/2 + vid[0]*Lx + kx*ix;
      py = -pi + vid[1]*Ly + ky*iy;
      IniCond();
      tred2(h, N, d, e);
      tqli(d, e, N, h);
      //	eigsrt(d, h, N);
      eigvec(h, V, Vt, N);
      sumN1 = 0;
      sumN2 = 0;
      sumM1 = 0;
      sumM2 = 0;
      for(k = 1; k <= N; k++)
	if(d[k] <= eF) { 
	  for(i = 1; i <= N; i++)
	  for(j = 1; j <= N; j++){
	    sumN1 += Vt[k][i]*N1[i][j]*V[j][k]; 
	    sumN2 += Vt[k][i]*N2[i][j]*V[j][k]; 
	    sumM1 += Vt[k][i]*M1[i][j]*V[j][k]; 
	    sumM2 += Vt[k][i]*M2[i][j]*V[j][k]; 
	  }
	}
      n1new += sumN1;
      n2new += sumN2;
      m1new += sumM1;
      m2new += sumM2;
    }
    MPI_Allreduce(&n1new,&n1newG,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&n2new,&n2newG,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&m1new,&m1newG,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&m2new,&m2newG,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    n1old = n1newG/NN;
    n2old = n2newG/NN;
    m1old = m1newG/NN;
    m2old = m2newG/NN;
    //    fprintf(p1,"%d  %f  %f  %f  %f  %f\n", sweep, n1new, n2new, m1new, m2new,eF);

  }

}
/*-----------------------------------------------------------
  Searching for Fermi surface.
-----------------------------------------------------------*/
void eSearch() {
  NE = 2*NN;
  Nh = N/2;
  NF = (int)(n*NE);

  hu = dmatrix(1,Nh,1,Nh);
  du = dvector(1,Nh);
  eu = dvector(1,Nh);
  eOrder = dvector(1,NE);

  count = 0;
  for(ix = -Nx/2; ix < Nx/2; ix++)
  for(iy = -Ny/2; iy < Ny/2; iy++) {
    px = kx*ix;
    py = ky*iy;
    IniCond();

    for(i = 1; i <= Nh; i++)
      for(j = 1; j <= Nh; j++) {
	hu[i][j] = h[i][j];
      }
      tred2(hu, Nh, du, eu);
      tqli(du, eu, Nh, hu);
      for(i = 1; i <= Nh; i++)
	eOrder[i + 4*count] = du[i];
      count++ ;
  }
  shell(NE, eOrder);
  eF = eOrder[NF];
}
/*------------------------------------------------------------
  Initial Conditions
------------------------------------------------------------*/
void IniCond() {
  double sign;

  n1 = n1old;
  n2 = n2old;
  m1 = m1old;
  m2 = m2old;
  
  ek11 = -2*t1*cos(px);
  ek12 = -2*t1*cos(py);
  ek21 = -2*t2*cos(py);
  ek22 = -2*t2*cos(px);
  ek3 =-4*t3*cos(px)*cos(py);
  ek4 = 4*t4*sin(px)*sin(py);

  ek1 = ek11 + ek21 + ek3 + U*n1 + 2*Uh*n2;
  ek2 = ek12 + ek22 + ek3 + U*n2 + 2*Uh*n1;

  ekq1 =-ek11 + ek21 - ek3 + U*n1 + 2*Uh*n2;
  ekq2 = ek12 - ek22 - ek3 + U*n2 + 2*Uh*n1;

  i1s = -(U*m1 + J*m2)/2;
  i2s = -(U*m2 + J*m1)/2;

  for(i = 1; i <= N; i++)
    for(j = 1; j <= N; j++){
      N1[i][j] = 0;
      N2[i][j] = 0;
      M1[i][j] = 0;
      M2[i][j] = 0;
      h[i][j]  = 0;
    }

  // non-zero elements of Hamiltonian matrix
  for(i = 0; i <=1; i++){
    sign = pow(-1,i);
    h[1+4*i][1+4*i] = ek1;
    h[1+4*i][2+4*i] = ek4;
    h[1+4*i][3+4*i] = sign*i1s;
    h[2+4*i][1+4*i] = h[1+4*i][2+4*i];
    h[2+4*i][2+4*i] = ek2;
    h[2+4*i][4+4*i] = sign*i2s;
    h[3+4*i][1+4*i] = h[1+4*i][3+4*i];
    h[3+4*i][3+4*i] = ekq1;
    h[3+4*i][4+4*i] = -ek4;
    h[4+4*i][2+4*i] = h[2+4*i][4+4*i];
    h[4+4*i][3+4*i] = h[3+4*i][4+4*i];
    h[4+4*i][4+4*i] = ekq2;
  }
    
  N1[1][1] = 0.5;
  N1[3][3] = 0.5;
  N1[5][5] = 0.5;
  N1[7][7] = 0.5;

  N2[2][2] = 0.5;
  N2[4][4] = 0.5;
  N2[6][6] = 0.5;
  N2[8][8] = 0.5;

  M1[1][3] = 1;
  M1[3][1] = 1;
  M1[5][7] = -1;
  M1[7][5] = -1;

  M2[2][4] = 1;
  M2[4][2] = 1;
  M2[6][8] = -1;
  M2[8][6] = -1;

  return;
}
