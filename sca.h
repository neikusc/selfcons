/*------------------------------------------------------------------
  Define parameters by self-consistent solving
------------------------------------------------------------------*/
#define pi 3.1415926
#define Nx 32
#define Ny 32
#define N 8
/*--------------------- parallel parameters ---------------------*/
int myid, nproc;
int vid[2];
int vproc[2] = {2,4}, nproc = 8;

int NN;
double Lx, Ly;    // horizontal and vertical local length 
double kx, ky;

/*------------------------------------------------------------------
  Define prototype functions
  1) Contruct Hamiltonian matrix and matrices define occupancy 
  numbers and staggered magnetization.
  2) Diagonalize symmetric matrix
------------------------------------------------------------------*/
void IniCond();
void eSearch();
void CompPara();
void tred2(double **, int , double *, double *); 
void tqli(double *, double *, int, double **);

/*------------------------------------------------------------------
  Define variables and matrices
------------------------------------------------------------------*/
int i, j ,k, count, sweep, Nmax;
int ix, iy;
int NE, NF;                              // NE: number of energy levels, NF: position of Fermi level
int Nh;                                         
double px, py, rBZ; 
double ek11, ek12, ek21, ek22, ek3, ek4; // interaction hopping coefficient
double ek1, ek2, ekq1, ekq2;
double i1s, i2s;                         // spin interaction
double sumN1, sumN2, sumM1, sumM2;
double n;
double n1, n2, m1, m2;                   // occupancy number and staggered number
double n1old, n2old, m1old, m2old;
double n1new, n2new, m1new, m2new;
double n1newG, n2newG, m1newG, m2newG;
double eF;                               // Fermi 
double U, J, Uh;

double **h, **hu;                        // Hamiltonian Matrix
double **N1, **N2, **M1, **M2;           // Matrices define occupancy numbers and staggered magnetization
double **V, **Vt;                        // EigenVectors define transformation
double *d, *e, *du, *eu;                 // EigenValues
double *eOrder;                          // Sort energy levels
