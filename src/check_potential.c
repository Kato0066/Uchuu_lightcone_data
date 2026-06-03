#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#define PI (3.141592653589793)
#define GRAVITY_G (1.0)
#define RHS_TINY (1.0e-30)

#define INDEX_2D(ix,iy,ngrid) ((iy) + (int64_t)(ngrid) * (ix))//1次元配列に直す

void read_delta_slice(const char *filename, double *delta, int ngrid)
{
  FILE *fp;
  double x, y, value;
  int64_t count = 0;
  int64_t datasize = (int64_t)ngrid * ngrid;

  fp = fopen(filename, "r");

  while (count < datasize &&
         fscanf(fp, "%lf %lf %lf", &x, &y, &value) == 3) {
    delta[count] = value;
    count++;
  }

  fclose(fp);
}

void read_potential_slices_binary(const char *filename,
                                  double *phi_m, double *phi_0, double *phi_p,
                                  int ngrid, int iz_slice)
{
  FILE *fp;
  int izp = (iz_slice + 1) % ngrid;
  int izm = (iz_slice - 1 + ngrid) % ngrid;//周期的境界条件

  fp = fopen(filename, "rb");

  for (int ix = 0; ix < ngrid; ix++) {
    for (int iy = 0; iy < ngrid; iy++) {
      int64_t im2 = INDEX_2D(ix, iy, ngrid);
      int64_t base = (int64_t)ngrid * ((int64_t)iy + (int64_t)ngrid * ix);

      fseek(fp, (base + izm) * (int64_t)sizeof(double), SEEK_SET);
      fread(&phi_m[im2], sizeof(double), 1, fp);

      fseek(fp, (base + iz_slice) * (int64_t)sizeof(double), SEEK_SET);
      fread(&phi_0[im2], sizeof(double), 1, fp);

      fseek(fp, (base + izp) * (int64_t)sizeof(double), SEEK_SET);
      fread(&phi_p[im2], sizeof(double), 1, fp);
    }
  }

  fclose(fp);
}

int main(int argc, char **argv)
{
  double *delta;
  double *phi_m;
  double *phi_0;
  double *phi_p;
  FILE *fout;
  int iz_slice;
  int ngrid;
  double boxsize;
  double scale_factor;
  double dx;
  int64_t slice_size;
  double max_rel_res = 0.0;

  if (argc != 8) {
    fprintf(stderr,
            "Usage: %s delta_slice.dat potential.bin output.dat iz_slice ngrid boxsize scale_factor\n",
            argv[0]);
    return 1;
  }

  iz_slice = atoi(argv[4]);
  ngrid = atoi(argv[5]);
  boxsize = atof(argv[6]);
  scale_factor = atof(argv[7]);
  dx = 1 / (double)ngrid;
  slice_size = (int64_t)ngrid * ngrid;

  delta = (double *)malloc(sizeof(double) * slice_size);
  phi_m = (double *)malloc(sizeof(double) * slice_size);
  phi_0 = (double *)malloc(sizeof(double) * slice_size);
  phi_p = (double *)malloc(sizeof(double) * slice_size);

  read_delta_slice(argv[1], delta, ngrid);
  read_potential_slices_binary(argv[2], phi_m, phi_0, phi_p,
                               ngrid, iz_slice);

  fout = fopen(argv[3], "w");

  for (int ix = 0; ix < ngrid; ix++) {
    int ixp = (ix + 1) % ngrid;//周期境界込み
    int ixm = (ix - 1 + ngrid) % ngrid;

    for (int iy = 0; iy < ngrid; iy++) {
      int iyp = (iy + 1) % ngrid;
      int iym = (iy - 1 + ngrid) % ngrid;

      int64_t im2 = INDEX_2D(ix, iy, ngrid);
      double x = dx * ((double)ix + 0.5);//セルの中心
      double y = dx * ((double)iy + 0.5);
      double lap_x;
      double lap_y;
      double lap_z;
      double laplacian;
      double rhs;
      double residual;
      double rel_residual;

      lap_x = (phi_0[INDEX_2D(ixp, iy, ngrid)]
               - 2.0 * phi_0[im2]
               + phi_0[INDEX_2D(ixm, iy, ngrid)]) / (dx * dx);
      lap_y = (phi_0[INDEX_2D(ix, iyp, ngrid)]
               - 2.0 * phi_0[im2]
               + phi_0[INDEX_2D(ix, iym, ngrid)]) / (dx * dx);
      lap_z = (phi_p[im2]
               - 2.0 * phi_0[im2]
               + phi_m[im2]) / (dx * dx);

      laplacian = lap_x + lap_y + lap_z;

      rhs = 4.0 * PI * scale_factor * scale_factor * delta[im2];
      residual = laplacian - rhs;

      if (fabs(rhs) > RHS_TINY) {
        rel_residual = fabs(residual) / fabs(rhs);
      } else {
        rel_residual = fabs(residual);
      }
      //0割回避

      fprintf(fout,
              "%14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e\n",
              x, y, delta[im2], phi_0[im2],
              laplacian, rhs, residual, rel_residual);
    }
    fprintf(fout, "\n");
  }

  fclose(fout);

  printf("# z-slice iz = %d\n", iz_slice);
  printf("# wrote %s\n", argv[3]);

  free(delta);
  free(phi_m);
  free(phi_0);
  free(phi_p);
  return 0;
}

#undef INDEX_2D
