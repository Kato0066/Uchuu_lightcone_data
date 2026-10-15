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

#define NOBS (8)

static const double OBS_POS[NOBS][3] = {
  {0.60,   0.50,   0.50},
  {0.70,   0.50,   0.50},
  {0.80,   0.50,   0.50},
  {0.90,   0.50,   0.50},
  {0.50,   0.70,   0.50},
  {0.50,   0.50,   0.70},
  {0.6414, 0.6414, 0.50},
  {0.6155, 0.6155, 0.6155}
};

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

void phi_const(double phi[N][N][N])
{
  for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){
      for(int k=0;k<N;k++){
        phi[i][j][k] = 0.00;
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
  const char *output_dir = "phi_const_001";
  char observer_filename[512];

  int iobs[NOBS];
  int jobs[NOBS];
  int kobs[NOBS];
  double xobs[NOBS];
  double yobs[NOBS];
  double zobs[NOBS];
  double robs[NOBS];
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

  if(argc == 1) {
    phi_const(phi);
  } else if(argc == 2) {
    output_dir = argv[1];
    phi_const(phi);
  } else if(argc == 4) {
    output_dir = argv[1];
    read_potential_binary(argv[2], phi);
    rho_bar = read_rho_bar(argv[3]);
    beta = calc_beta(rho_bar);
  }
  snprintf(observer_filename, sizeof(observer_filename),
           "%s/fixed_point.dat", output_dir);

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

  for(int n=0; n<NOBS; n++) {
    iobs[n] = (int)(OBS_POS[n][0] / dx);
    jobs[n] = (int)(OBS_POS[n][1] / dx);
    kobs[n] = (int)(OBS_POS[n][2] / dx);


    xobs[n] = dx*((double)iobs[n] + 0.5);
    yobs[n] = dx*((double)jobs[n] + 0.5);
    zobs[n] = dx*((double)kobs[n] + 0.5);

    robs[n] = sqrt(
                   (xobs[n]-0.5)*(xobs[n]-0.5)
                   + (yobs[n]-0.5)*(yobs[n]-0.5)
                   + (zobs[n]-0.5)*(zobs[n]-0.5));
  }

  fp_observer = fopen("fixed_point.dat", "w");

  fprintf(fp_observer,
          "# 1:tau 2:id 3:x 4:y 5:z 6:r "
          "7:F0 8:D 9:normalized_phi 10:F=F0+D\n");

  for(int n=0; n<NOBS; n++) {
    fprintf(fp_observer,
            "# observer %d: index=(%d,%d,%d), "
            "position=(%.8e,%.8e,%.8e), r=%.8e\n",
            n, iobs[n], jobs[n], kobs[n],
            xobs[n], yobs[n], zobs[n], robs[n]);

    printf("# observer %d: index=(%d,%d,%d), "
           "position=(%.8e,%.8e,%.8e), r=%.8e\n",
           n, iobs[n], jobs[n], kobs[n],
           xobs[n], yobs[n], zobs[n], robs[n]);
  }

  for(int n=0; n<NOBS; n++) {
    fprintf(fp_observer,
            "%14.6e %d "
            "%14.6e %14.6e %14.6e %14.6e "
            "%14.6e %14.6e %14.6e %14.6e\n",
            tnow, n,
            xobs[n], yobs[n], zobs[n], robs[n],
            F01[iobs[n]][jobs[n]][kobs[n]],
            D1[iobs[n]][jobs[n]][kobs[n]],
            beta*phi[iobs[n]][jobs[n]][kobs[n]],
            F01[iobs[n]][jobs[n]][kobs[n]]
            + D1[iobs[n]][jobs[n]][kobs[n]]);
  }

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

    for(int n=0; n<NOBS; n++) {
      fprintf(fp_observer,
              "%14.6e %d "
              "%14.6e %14.6e %14.6e %14.6e "
              "%14.6e %14.6e %14.6e %14.6e\n",
              tnow, n,
              xobs[n], yobs[n], zobs[n], robs[n],
              F01[iobs[n]][jobs[n]][kobs[n]],
              D1[iobs[n]][jobs[n]][kobs[n]],
              beta*phi[iobs[n]][jobs[n]][kobs[n]],
              F01[iobs[n]][jobs[n]][kobs[n]]
              + D1[iobs[n]][jobs[n]][kobs[n]]);
    }

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
