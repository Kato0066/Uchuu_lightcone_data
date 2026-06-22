#include <stdio.h>
#include <math.h>

#define N (128)

int main(int argc, char **arg){
  int i, j, k;
  double nu_cfl = 0.5;
  double c = 1.0;
  double dx = 1.0 / (double)N;
  double dt = nu_cfl*dx/c;
  double PI = 3.141592653589793;

  double tend;
  int istep=0;

  static double U0[N][N][N];
  static double U1[N][N][N];
  static double U2[N][N][N];
  static double phi[N][N][N];
  static double A[N][N][N]

    //初期条件
    for(i=0;i<N;i++){
      for(j=0;j<N;j++){
	for(k=0;k<N;k++){
	  double x = (i+0.5)*dx;
	  double y = (j+0.5)*dx;
	  double z = (k+0.5)*dx;
	  U0[i][j]= exp(-100.0*((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)+(z-0.5)*(z-0.5)));
	  //U0[i][j]= sin(2.0*PI*x);
	}
      }
    }


  double phi[N][N][N], A[N][N][N];
  for(i=0; i<N; i++){
    for(j=0; j<N; j++){
      for(k=0; k<N; k++){
        double x = (i+0.5)*dx;
        double y = (j+0.5)*dx;
        double z = (k+0.5)*dx;
        phi[i][j][k] = 0.0;
        A[i][j][k] = 1.0/(1.0-4.0*phi[i][j][k])*(nu_cfl*nu_cfl);

        int ip = i+1 > N-1 ? i+1-N : i+1;
        int im = i-1 < 0 ? i-1+N : i-1;
        int jp = j+1 > N-1 ? j+1-N : j+1;
        int jm = j-1 < 0 ? j-1+N : j-1;
        int kp = k+1 > N-1 ? k+1-N : k+1;
        int km = k-1<0 ? k-1+N : k-1;
        U1[i][j][k] = U0[i][j][k] + 0.5 * A[i][j][k] * (U0[ip][j][k] + U0[im][j][k] + U0[i][jp][k] + U0[i][jm][k] + U0[i][j][kp] + U0[i][j][km] - 6.0 * U0[i][j][k]);
        //U1[i] = U0[i];
      }
    }
  }



  double tnow = 0.0;
  while(tnow < tend) {
    for(i=0;i<N;i++){
      for(j=0;j<N;j++){
        for(k=0;k<N;k++){
          tnow += dt;
          istep += 1;
          times[istep] = tnow;
          int ip = i+1 > N-1 ? i+1-N : i+1;
          int im = i-1 < 0 ? i-1+N : i-1;
          int jp = j+1 > N-1 ? j+1-N : j+1;
          int jm = j-1 < 0 ? j-1+N : j-1;
          int kp = k+1 > N-1 ? k+1-N : k+1;
          int km = k-1<0 ? k-1+N : k-1;
          U2[i][j][k] = 2.0 * U1[i][j][k] - U0[i][j][k] + A[i][j][k] * (U1[ip][j][k] + U1[im][j][k] + U1[i][jp][k] + U1[i][jm][k] + U1[i][j][kp] + U1[i][j][km] - 6.0 * U1[i][j][k]);

          U0[i][j][k] = U1[i][j][k];
          U1[i][j][k] = U2[i][j][k];
        }
      }
    }
  }
