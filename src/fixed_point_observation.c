#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N (512)
#define PI (3.141592653589793)
#define BOXSIZE (2000.0)
#define GRAVITY_G (4.30091727003628e-9)
#define SPEED_OF_LIGHT (299792.458)
#define GADGET_UNIT_MASS_IN_MSUN (1.0e10)
#define WAVE_AMPLITUDE (1.0)

#define OBS_X (0.75)
#define OBS_Y (0.5)
#define OBS_Z (0.5)


double source_Q(int i, int j, int k, double t,
                double c, double dx)
{
  int ic1 = N/2 -1;
  int ic2 = N/2;
  double wave_number = 100 * PI;
  double omega = c * wave_number;
  double T = WAVE_AMPLITUDE * cos(omega * t);

  if((i == ic1 || i == ic2)
     && (j == ic1 || j == ic2)
     && (k == ic1 || k ==ic2)){
    return T / (8 * dx * dx * dx);
  }

  return 0.0;
}

void phi_zero(double phi[N][N][N])
{
  for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){
      for(int k=0;k<N;k++){
	phi[i][j][k] = 0.0;
      }
    }
  }
}

void read_potential_binary(const char *filename, double phi[N][N][N])
{
  FILE *fp;

  fp = fopen(filename,"rb");
  fread(&phi[0][0][0],sizeof(double), N*N*N,fp);
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
  double tend = 1.0;
  int istep = 1;
  double beta = 1.0;
  double rho_bar = 0.0;

  /* 指定した座標に最も近いセルを観測点にする。 */
  int iobs = (int)(OBS_X / dx);
  int jobs = (int)(OBS_Y / dx);
  int kobs = (int)(OBS_Z / dx);
  double xobs;
  double yobs;
  double zobs;
  FILE *fp_observer;

  double (*D0)[N][N];
  double (*D1)[N][N];
  double (*D2)[N][N];
  double (*phi)[N][N];
  double (*F00)[N][N];
  double (*F01)[N][N];
  double (*F02)[N][N];

  D0 = malloc(sizeof(double) * N * N * N);
  D1 = malloc(sizeof(double) * N * N * N);
  D2 = malloc(sizeof(double) * N * N * N);
  phi = malloc(sizeof(double) * N * N * N);
  F00 = malloc(sizeof(double) * N * N * N);
  F01 = malloc(sizeof(double) * N * N * N);
  F02 = malloc(sizeof(double) * N * N * N);

  if(argc == 3) {
    read_potential_binary(argv[1], phi);
    rho_bar = read_rho_bar(argv[2]);
    beta = calc_beta(rho_bar);
  } else {
    phi_zero(phi);
  }

  printf("# rho_bar = %.8e\n", rho_bar);
  printf("# beta = %.8e\n", beta);
  printf("# normalized phi at source = %.8e\n",
         beta*phi[N/2][N/2][N/2]);

#pragma omp parallel for collapse(3) schedule(auto) private(i,j,k)
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      for(k=0;k<N;k++){
        D0[i][j][k] = 0.0;
        D1[i][j][k] = 0.0;
        F00[i][j][k] = 0.0;
        F01[i][j][k] = 0.0;
      }
    }
  }

  if(iobs < 0 || iobs >= N
     || jobs < 0 || jobs >= N
     || kobs < 0 || kobs >= N) {
    fprintf(stderr,
            "OBS_X, OBS_Y and OBS_Z must satisfy 0 <= OBS_* < 1\n");
    return 1;
  }

  xobs = dx*((double)iobs + 0.5);
  yobs = dx*((double)jobs + 0.5);
  zobs = dx*((double)kobs + 0.5);

  fp_observer = fopen("fixed_point.dat", "w");
  if(fp_observer == NULL) {
    fprintf(stderr, "cannot open fixed-point output file\n");
    return 1;
  }

  fprintf(fp_observer, "# x = %.8e\n", xobs);
  fprintf(fp_observer, "# y = %.8e\n", yobs);
  fprintf(fp_observer, "# z = %.8e\n", zobs);
  fprintf(fp_observer, "# normalized_phi = %.8e\n",
          beta*phi[iobs][jobs][kobs]);
  fprintf(fp_observer,
          "# 1:tau 2:x 3:y 4:z 5:F0 6:D 7:normalized_phi\n");
  fprintf(fp_observer,
          "%14.6e %14.6e %14.6e %14.6e "
          "%14.6e %14.6e %14.6e\n",
          tnow, xobs, yobs, zobs,
          F01[iobs][jobs][kobs],
          D1[iobs][jobs][kobs],
          beta*phi[iobs][jobs][kobs]);

  printf("# fixed observer: index=(%d,%d,%d), position=(%.8e,%.8e,%.8e)\n",
         iobs, jobs, kobs, xobs, yobs, zobs);

  while(tnow < tend) {
#pragma omp parallel for collapse(3) schedule(auto) private(i,j,k)
    for(i=0;i<N;i++){
      for(j=0;j<N;j++){
        for(k=0;k<N;k++){
          double normalized_phi = beta*phi[i][j][k];
          double denominator = 1.0 - 4.0 * normalized_phi;

          int ip = i+1 > N-1 ? i+1-N : i+1;
          int im = i-1 < 0 ? i-1+N : i-1;
          int jp = j+1 > N-1 ? j+1-N : j+1;
          int jm = j-1 < 0 ? j-1+N : j-1;
          int kp = k+1 > N-1 ? k+1-N : k+1;
          int km = k-1 < 0 ? k-1+N : k-1;

          double lap_F0;
          double Q;
          double S;

          lap_F0 =
            (F01[ip][j][k] + F01[im][j][k]
             + F01[i][jp][k] + F01[i][jm][k]
             + F01[i][j][kp] + F01[i][j][km]
             - 6.0 * F01[i][j][k]) / (dx*dx);

          Q = source_Q(i, j, k, tnow, c, dx);
          S = lap_F0 + Q;

          F02[i][j][k] = 2.0 * F01[i][j][k] - F00[i][j][k]
            + dt*dt*S;

          D2[i][j][k] = 2.0 * D1[i][j][k] - D0[i][j][k]
            + dt*dt/(denominator*dx*dx)
            * (D1[ip][j][k] + D1[im][j][k]
               + D1[i][jp][k] + D1[i][jm][k]
               + D1[i][j][kp] + D1[i][j][km]
               - 6.0 * D1[i][j][k])
            + dt*dt/denominator
            * 4.0*normalized_phi*S;
        }
      }
    }

#pragma omp parallel for collapse(3) schedule(auto) private(i,j,k)
    for(i=0;i<N;i++){
      for(j=0;j<N;j++){
        for(k=0;k<N;k++){
          D0[i][j][k] = D1[i][j][k];
          D1[i][j][k] = D2[i][j][k];
          F00[i][j][k] = F01[i][j][k];
          F01[i][j][k] = F02[i][j][k];
        }
      }
    }

    tnow += dt;
    istep += 1;

    /* 観測点の時間発展を毎ステップ記録する。 */
    fprintf(fp_observer,
            "%14.6e %14.6e %14.6e %14.6e "
            "%14.6e %14.6e %14.6e\n",
            tnow, xobs, yobs, zobs,
            F01[iobs][jobs][kobs],
            D1[iobs][jobs][kobs],
            beta*phi[iobs][jobs][kobs]);

    if(istep%100==0 || tnow>=tend) {
      printf("# step %d: tau = %.8e\n", istep, tnow);
      fflush(stdout);
      fflush(fp_observer);
    }

  }

  fclose(fp_observer);

  free(D0);
  free(D1);
  free(D2);
  free(phi);
  free(F00);
  free(F01);
  free(F02);

  return 0;
}
