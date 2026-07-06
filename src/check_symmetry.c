#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N (512)
#define NSTEP (200)

void read_wave_binary(const char *filename, double u[N][N][N])
{
  FILE *fp;

  fp = fopen(filename, "rb");
  fread(&u[0][0][0], sizeof(double), N*N*N, fp);
  fclose(fp);
}

int main(int argc, char **argv)
{
  double (*u)[N][N];
  FILE *fp;
  double dx = 1.0 / (double)N;
  int ic = N / 2;
  int jc = N / 2;
  int kc = N / 2;

  u = malloc(sizeof(double) * N * N * N);

  read_wave_binary(argv[1], u);

  fp = fopen(argv[2], "w");
  fprintf(fp, "# direction r u\n");

  for(int m=0;m<NSTEP;m++) {
    fprintf(fp, "x %14.6e %14.6e\n", dx * (double)m, u[ic + m][jc][kc]);
}
fprintf(fp, "\n");

for(int m=0;m<NSTEP;m++) {
  fprintf(fp, "y %14.6e %14.6e\n", dx * (double)m, u[ic][jc + m][kc]);
 }
fprintf(fp, "\n");

for(int m=0;m<NSTEP;m++) {
  fprintf(fp, "z %14.6e %14.6e\n", dx * (double)m, u[ic][jc][kc + m]);
 }
fprintf(fp, "\n");

for(int m=0;m<NSTEP;m++) {
  fprintf(fp, "xy %14.6e %14.6e\n", sqrt(2.0) * dx * (double)m, u[ic + m][jc + m][kc]);
 }
fprintf(fp, "\n");

for(int m=0;m<NSTEP;m++) {
  fprintf(fp, "xz %14.6e %14.6e\n",sqrt(2.0) * dx * (double)m, u[ic + m][jc][kc + m]);
 }
fprintf(fp, "\n");

for(int m=0;m<NSTEP;m++) {
  fprintf(fp, "yz %14.6e %14.6e\n",sqrt(2.0) * dx * (double)m, u[ic][jc + m][kc + m]);
 }
fprintf(fp, "\n");

for(int m=0;m<NSTEP;m++) {
  fprintf(fp, "xyz %14.6e %14.6e\n", sqrt(3.0) * dx * (double)m, u[ic + m][jc + m][kc + m]);
 }
fprintf(fp, "\n");

fclose(fp);
free(u);

return 0;
}
