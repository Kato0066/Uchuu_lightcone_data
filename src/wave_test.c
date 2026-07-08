#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N (512)
#define PI (3.141592653589793)
#define BOXSIZE (2000.0)
#define GRAVITY_G (4.30091727003628e-9)
#define SPEED_OF_LIGHT (299792.458)
#define GADGET_UNIT_MASS_IN_MSUN (1.0e10)

double initial_wave(double x, double y, double z, double t, double c, double sigma)
{
  double xc = 0.5;
  double yc = 0.5;
  double zc = 0.5;
  double r0 = 0.05;

  double r = sqrt((x-xc)*(x-xc) + (y-yc)*(y-yc) + (z-zc)*(z-zc));

  if(r < 1.0e-12) return 0.0;

  return (r0/r) * exp(-0.5 * (r - r0 - c*t) * (r - r0 - c*t) / (sigma * sigma));


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

void output_wave_binary(const char *filename, double u[N][N][N])
{
  FILE *fp;

  fp = fopen(filename, "wb");
  fwrite(&u[0][0][0], sizeof(double), N*N*N, fp);
  fclose(fp);
}

double read_rho_bar(const char *filename)
{
  FILE *fp;
  double rho_bar;

  fp = fopen(filename, "r");
  fscanf(fp, "%lf", &rho_bar);
  fclose(fp);

  return rho_bar;
}

double calc_beta(double rho_bar)
{
  double total_mass_msun;

  total_mass_msun = rho_bar * (double)N * (double)N * (double)N * GADGET_UNIT_MASS_IN_MSUN;

  return GRAVITY_G * total_mass_msun
    / (BOXSIZE * SPEED_OF_LIGHT * SPEED_OF_LIGHT);
}


int main(int argc, char **argv){
  int i, j, k;
  double nu_cfl = 0.25;
  double c = 1.0;
  double dx = 1.0 / (double)N;
  double dt = nu_cfl*dx/c;
  //double PI = 3.141592653589793;
  double tnow = 0.0;
  double tend = 1.0;
  int istep = 0;
  int jmid = N / 2;
  int kmid = N / 2;
  double sigma = 5*dx;
  double beta = 1.0;
  double rho_bar = 0.0;

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

  if(argc == 3) {
    read_potential_binary(argv[1], phi);
    rho_bar = read_rho_bar(argv[2]);
    beta = calc_beta(rho_bar);
  } else {
    phi_zero(phi);
  }

#pragma omp parallel for collapse(3) schedule(auto) private(i,j,k)
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      for(k=0;k<N;k++){
        double x = (i+0.5)*dx;
        double y = (j+0.5)*dx;
        double z = (k+0.5)*dx;

        U0[i][j][k] = initial_wave(x, y, z, 0.0, c, sigma);
        U1[i][j][k] = initial_wave(x, y, z, dt, c, sigma);
      }
    }
  }

#if 0
#pragma omp parallel for collapse(3) schedule(auto) private(i,j,k)
  for(i=0; i<N; i++){
    for(j=0; j<N; j++){
      for(k=0; k<N; k++){
        A[i][j][k] = 1.0/(1.0-4.0*beta*phi[i][j][k])*(nu_cfl*nu_cfl);

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


  printf("tendの値を入力してください。\n");
  scanf("%lf", &tend);
#endif


  FILE *output_wave_slice;
  output_wave_slice = fopen("wave_slices_cfl025/wave_slice_final.dat", "w");

  while(tnow < tend) {
#pragma omp parallel for collapse(3) schedule(auto) private(i,j,k)
    for(i=0;i<N;i++){
      for(j=0;j<N;j++){
        for(k=0;k<N;k++){
          A[i][j][k] = 1.0/(1.0-4.0*beta*phi[i][j][k])*(nu_cfl*nu_cfl);

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

#pragma omp parallel for collapse(3) schedule(auto) private(i,j,k)
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

    if(istep%128==0) {
      char binary_name[120];

      sprintf(binary_name, "wave_binary_cfl025/wave_step%04d_t%06.4f.bin", istep, tnow);
      output_wave_binary(binary_name, U2);
    }

    if(istep%128==0) {
      FILE *fp;
      FILE *fp_line;
      char name[120];
      char line_name[120];

      sprintf(name, "wave_slices_cfl025/wave_slice-%04d_t%06.4f.dat", istep, tnow);
      fp = fopen(name, "w");
      for(i=0;i<N;i++) {
        for(j=0;j<N;j++) {
          fprintf(fp, "%12.4e %12.4e %12.4e\n",
                  dx*(double)i, dx*(double)j, U2[i][j][kmid]);
        }
        fprintf(fp, "\n");
      }
      fclose(fp);

      sprintf(line_name, "wave_slices_cfl025/wave_line-%04d_t%06.4f.dat", istep, tnow);
      fp_line = fopen(line_name, "w");
      for(i=0;i<N;i++) {
        fprintf(fp_line, "%12.4e %12.4e\n",
                dx*((double)i+0.5), U1[i][jmid][kmid]);
      }
      fclose(fp_line);
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

  {
    FILE *output_wave_line;
    char final_line_name[128];

    sprintf(final_line_name, "wave_slices_cfl025/wave_line_final.dat");
    output_wave_line = fopen(final_line_name, "w");

    for(i=0;i<N;i++) {
      fprintf(output_wave_line, "%12.4e %12.4e\n",
              dx*((double)i+0.5), U1[i][jmid][kmid]);
    }

    fclose(output_wave_line);
  }

  free(U0);
  free(U1);
  free(U2);
  free(phi);
  free(A);

  return 0;

}
