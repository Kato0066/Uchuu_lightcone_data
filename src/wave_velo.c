#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N (512)

void read_wave_binary(const char *filename, double u[N][N][N])
{
  FILE *fp;

  fp = fopen(filename, "rb");
  fread(&u[0][0][0], sizeof(double), N*N*N, fp);
  fclose(fp);
}

int main(void)
{
  double (*u)[N][N];
  FILE *fp_out;
  double nu_cfl = 0.5;
  double c = 1.0;
  double dx = 1.0 / (double)N;
  double dt = nu_cfl * dx / c;
  double x_prev = 0.0;
  double t_prev = 0.0;
  int jmid = N / 2;
  int kmid = N / 2;

  u = malloc(sizeof(double) * N * N * N);
  fp_out = fopen("peak_velocity.dat", "w");

  fprintf(fp_out, "# step t x_peak U_peak velocity\n");

  for(int step=128;step<=1024;step+=128) {
    char filename[64];
    double max_u = -INFINITY;
    double x_peak = 0.0;
    double velocity = 0.0;
    double t = dt * (double)step;

    sprintf(filename, "wave_binary/wave_step%05d.bin", step);
    read_wave_binary(filename, u);

    for(int i=N/2;i<N;i++) {
      if(u[i][jmid][kmid] > max_u) {
        max_u = u[i][jmid][kmid];
        x_peak = ((double)i + 0.5) * dx;
      }
    }

    if(step > 128) {
      velocity = (x_peak - x_prev) / (t - t_prev);
    }

    fprintf(fp_out, "%d %14.6e %14.6e %14.6e %14.6e\n",
            step, t, x_peak, max_u, velocity);

    x_prev = x_peak;
    t_prev = t;
  }

  fclose(fp_out);
  free(u);

  return 0;
}
