#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N (128)

double initial_wave(double x, double y, double z)
{
  return exp(-100.0*((x-0.5)*(x-0.5)
                     + (y-0.5)*(y-0.5)
                     + (z-0.5)*(z-0.5)));
}

void phi_zero(double phi[N][N][N])
{
  for(int i=0;i<N;i++) {
    for(int j=0;j<N;j++) {
      for(int k=0;k<N;k++) {
        phi[i][j][k] = 0.0;
      }
    }
  }
}

void read_potential_binary(const char *filename, double phi[N][N][N])
{
  FILE *fp;

  fp = fopen(filename, "rb");
  fread(&phi[0][0][0], sizeof(double), N*N*N, fp);
  fclose(fp);
}



int main(int argc, char **arg){
  int i, j, k;
  double nu_cfl = 0.5;
  double c = 1.0;
  double dx = 1.0 / (double)N;
  double dt = nu_cfl*dx/c;
  double PI = 3.141592653589793;
  double tnow = 0.0;
  double tend = 0.7;
  int istep = 0;
  int kmid = N / 2;

  double (*U0)[N][N];
  double (*U1)[N][N];
  double (*U2)[N][N];
  double (*phi)[N][N];
  double (*A)[N][N];

  U0 = malloc(sizeof(double) * N * N * N);
  U1 = malloc(sizeof(double) * N * N * N);
  U2 = malloc(sizeof(double) * N * N * N);
  phi = malloc(sizeof(double) * N * N * N);
  A = malloc(sizeof(double) * N * N * N);

  if(argc == 2) {
    read_potential_binary(argv[1], phi);
  } else {
    phi_zero(phi);
  }

  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      for(k=0;k<N;k++){
        double x = (i+0.5)*dx;
        double y = (j+0.5)*dx;
        double z = (k+0.5)*dx;

        U0[i][j][k] = initial_wave(x, y, z);
      }
    }
  }

  for(i=0; i<N; i++){
    for(j=0; j<N; j++){
      for(k=0; k<N; k++){
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

#if 0
  printf("tendの値を入力してください。\n");
  scanf("%lf", &tend);
#endif


  FILE *output_wave_slice;
  output_wave_slice = fopen("wave_slices/wave_slice_final.dat", "w");

  while(tnow < tend) {
    for(i=0;i<N;i++){
      for(j=0;j<N;j++){
        for(k=0;k<N;k++){
          int ip = i+1 > N-1 ? i+1-N : i+1;
          int im = i-1 < 0 ? i-1+N : i-1;
          int jp = j+1 > N-1 ? j+1-N : j+1;
          int jm = j-1 < 0 ? j-1+N : j-1;
          int kp = k+1 > N-1 ? k+1-N : k+1;
          int km = k-1 < 0 ? k-1+N : k-1;

          U2[i][j][k] =2.0 * U1[i][j][k] - U0[i][j][k] + A[i][j][k] * (U1[ip][j][k] + U1[im][j][k] + U1[i][jp][k] + U1[i][jm][k] + U1[i][j][kp] + U1[i][j][km] - 6.0 * U1[i][j][k]);
        }
      }
    }

    for(i=0;i<N;i++){
      for(j=0;j<N;j++){
        for(k=0;k<N;k++){
          U0[i][j][k] = U1[i][j][k];
          U1[i][j][k] = U2[i][j][k];
        }
      }
    }

    tnow += dt;
    istep += 1;

    if(istep%5==0) {
      FILE *fp;
      char name[120];

      sprintf(name, "wave_slices/wave_slice-%03d.dat", istep);
      fp = fopen(name, "w");
      for(i=0;i<N;i++) {
        for(j=0;j<N;j++) {
          fprintf(fp, "%12.4e %12.4e %12.4e\n",
                  dx*(double)i, dx*(double)j, U2[i][j][kmid]);
        }
        fprintf(fp, "\n");
      }
      fclose(fp);
    }
  }

  for(i=0;i<N;i++) {
    for(j=0;j<N;j++) {
      fprintf(output_wave_slice, "%12.4e %12.4e %12.4e\n",
              dx*(double)i, dx*(double)j, U2[i][j][kmid]);
    }
    fprintf(output_wave_slice, "\n");
  }
  fclose(output_wave_slice);

  free(U0);
  free(U1);
  free(U2);
  free(phi);
  free(A);

  return 0;

}
