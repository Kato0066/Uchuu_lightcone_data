#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#define N (512)
#define PI (3.141592653589793)
#define BOXSIZE (2000.0)
#define GRAVITY_G (4.30091727003628e-9)
#define SPEED_OF_LIGHT (299792.458)
#define GADGET_UNIT_MASS_IN_MSUN (1.0e10)

#define WAVE_AMPLITUDE (1.0)


double source_E(double x, double y, double z, double t,
                double c, double normalized_phi)
{
  double xs = 0.5;
  double ys = 0.5;
  double zs = 0.5;
  double wave_number = 100.0 * PI;
  double omega = c * wave_number;

  double r = sqrt((x-xs)*(x-xs)
                  + (y-ys)*(y-ys)
                  + (z-zs)*(z-zs));

  if(r < 1.0e-12) return 0.0;

  return 4.0 * WAVE_AMPLITUDE * omega * omega * normalized_phi / r
    * cos(wave_number * r - omega * t);
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

void output_D_binary(const char *filename, double D[N][N][N])
{
  FILE *fp;

  fp = fopen(filename, "wb");
  fwrite(&D[0][0][0], sizeof(double complex), N*N*N, fp);
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

  total_mass_msun = rho_bar * (double)N * (double)N * (double)N
    * GADGET_UNIT_MASS_IN_MSUN;

  return GRAVITY_G * total_mass_msun
    / (BOXSIZE * SPEED_OF_LIGHT * SPEED_OF_LIGHT);
}

int main(int argc, char **argv){
  int i, j, k;
  double nu_cfl = 0.25;
  double c = 1.0;
  double dx = 1.0 / (double)N;
  double dt = nu_cfl*dx/c;
  double tnow = dt;
  double tend = 0.04;
  int istep = 1;
  int jmid = N / 2;
  int kmid = N / 2;
  double beta = 1.0;
  double rho_bar = 0.0;

  double (*D0)[N][N];
  double (*D1)[N][N];
  double (*D2)[N][N];
  double (*phi)[N][N];
  double (*B)[N][N];

  D0 = malloc(sizeof(double) * N * N * N);
  D1 = malloc(sizeof(double) * N * N * N);
  D2 = malloc(sizeof(double) * N * N * N);
  phi = malloc(sizeof(double) * N * N * N);
  B = malloc(sizeof(double) * N * N * N);

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
        D0[i][j][k] = 0.0;
        D1[i][j][k] = 0.0;
      }
    }
  }

  FILE *output_D_slice;
  output_D_slice = fopen("D_slices_cfl025/D_slice_final.dat", "w");

  while(tnow < tend) {
#pragma omp parallel for collapse(3) schedule(auto) private(i,j,k)
    for(i=0;i<N;i++){
      for(j=0;j<N;j++){
        for(k=0;k<N;k++){
          double x = (i+0.5)*dx;
          double y = (j+0.5)*dx;
          double z = (k+0.5)*dx;
          double normalized_phi = beta*phi[i][j][k];
          double E = source_E(x, y, z, tnow, c, normalized_phi);

          B[i][j][k] = 1.0/(1.0-4.0*normalized_phi)
            *(nu_cfl*nu_cfl);

          int ip = i+1 > N-1 ? i+1-N : i+1;
          int im = i-1 < 0 ? i-1+N : i-1;
          int jp = j+1 > N-1 ? j+1-N : j+1;
          int jm = j-1 < 0 ? j-1+N : j-1;
          int kp = k+1 > N-1 ? k+1-N : k+1;
          int km = k-1 < 0 ? k-1+N : k-1;

          D2[i][j][k] = 2.0 * D1[i][j][k] - D0[i][j][k]
            + B[i][j][k]
            * (D1[ip][j][k] + D1[im][j][k]
               + D1[i][jp][k] + D1[i][jm][k]
               + D1[i][j][kp] + D1[i][j][km]
               - 6.0 * D1[i][j][k])
            - dt*dt/(1.0-4.0*normalized_phi)*E;
        }
      }
    }

#pragma omp parallel for collapse(3) schedule(auto) private(i,j,k)
    for(i=0;i<N;i++){
      for(j=0;j<N;j++){
        for(k=0;k<N;k++){
          D0[i][j][k] = D1[i][j][k];
          D1[i][j][k] = D2[i][j][k];
        }
      }
    }

    tnow += dt;
    istep += 1;

    printf("# step %d: tau = %.8e\n",istep, tnow);
    fflush(stdout);

    if(istep%128==0) {
      char binary_name[120];

      sprintf(binary_name,
              "D_binary_cfl025/D_step%04d_t%06.4f.bin",
              istep, tnow);
      output_D_binary(binary_name, D2);
    }

    if(istep%128==0) {
      FILE *fp;
      FILE *fp_line;
      char name[120];
      char line_name[120];

      sprintf(name,
              "D_slices_cfl025/D_slice-%04d_t%06.4f.dat",
              istep, tnow);
      fp = fopen(name, "w");
      for(i=0;i<N;i++) {
        for(j=0;j<N;j++) {
          fprintf(fp, "%12.4e %12.4e %12.4e %12.4e\n",
                  dx*(double)i, dx*(double)j,
                  creal(D2[i][j][kmid]),
                  cimag(D2[i][j][kmid]));
        }
        fprintf(fp, "\n");
      }
      fclose(fp);

      sprintf(line_name,
              "D_slices_cfl025/D_line-%04d_t%06.4f.dat",
              istep, tnow);
      fp_line = fopen(line_name, "w");
      for(i=0;i<N;i++) {
        fprintf(fp_line, "%12.4e %12.4e %12.4e\n",
                dx*((double)i+0.5),
                creal(D1[i][jmid][kmid]),
                cimag(D1[i][jmid][kmid]));
      }
      fclose(fp_line);
    }
  }

  for(i=0;i<N;i++) {
    for(j=0;j<N;j++) {
      fprintf(output_D_slice, "%12.4e %12.4e %12.4e\n",
              dx*((double)i+0.5), dx*((double)j+0.5),
              D2[i][j][kmid]);
    }
    fprintf(output_D_slice, "\n");
  }
  fclose(output_D_slice);

  {
    FILE *output_D_line;
    char final_line_name[128];

    sprintf(final_line_name, "D_slices_cfl025/D_line_final.dat");
    output_D_line = fopen(final_line_name, "w");

    for(i=0;i<N;i++) {
      fprintf(output_D_line, "%12.4e %12.4e\n",
              dx*((double)i+0.5),
              D1[i][jmid][kmid]);
    }

    fclose(output_D_line);
  }

  free(D0);
  free(D1);
  free(D2);
  free(phi);
  free(B);

  return 0;
}
