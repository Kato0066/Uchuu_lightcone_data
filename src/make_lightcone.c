/*Copyright (c) 2020 Tomoaki Ishiyama
  MIT License*/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>
#include <fftw3.h>


const double Gadget_UnitLength_in_Mpc = 1.0;            // 1Mpc
const double Gadget_UnitMass_in_Msun = 1.0e10;          // 1e10 Msun
const double Gadget_UnitVelocity_in_cm_per_s = 1e5;     //  1 km/sec

#define IO_CACHE_SIZE (2097152)
#define NGRID (512)
#define RHO(a,b,c) (mesh_density[(c) + (ngrid+2) * ((b) + ngrid * (a))])
#define BOXSIZE (2000.0)
#define PI (3.141592653589793)
#define SQR(x) ((x)*(x)) //square2乗
#define CUBE(x) ((x)*(x)*(x))
#define QUART(x) ((x)*(x)*(x)*(x))
#define GRAVITY_G (1.0) //後で計算



typedef struct Particle{
  double mass;
  double r[3];
} Particle;

typedef struct GadgetHeader {
  int      Npart[6];
  double   Massarr[6];
  double   Time;
  double   Redshift;
  int      FlagSfr;
  int      FlagFeedback;
  unsigned int Nall[6];
  int      FlagCooling;
  int      NumFiles;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam;
  int      Flag_StellarAge;
  int      Flag_Metals;
  unsigned int NallHW[6];
  int      flag_entr_ics;
  char     unused[256 - 6*4 - 6*8 - 2*8 - 2*4 - 6*4 - 2*4 - 4*8 - 9*4];
} GadgetHeader;



int get_gadget_npart(const char *filename)
{
  FILE *fin;
  int blksize;
  GadgetHeader gadget_header;

  fin = fopen(filename, "rb");
  if (fin == NULL) {
    fprintf(stderr, "cannot open %s\n", filename);
    return -1;
  }

  fread( &blksize, 4, 1, fin);
  //  fprintf( stderr, "Reading Gadget file blksize= %d\n", blksize);
  fread( &gadget_header, sizeof(gadget_header), 1, fin);
  fread( &blksize, 4, 1, fin);
  //  fprintf( stderr, "End Reading Gadget file blksize= %d\n", blksize);

  fclose(fin);

  return gadget_header.Npart[1];
}

int read_gadget_ptcl(const char *filename, Particle *ptcl)
{
  FILE *fin;
  int blksize;
  int npart;
  int i;
  GadgetHeader gadget_header;
  float *cache;

  int32_t nmemory = IO_CACHE_SIZE;

  fin = fopen(filename, "rb");
  if (fin == NULL) {
    fprintf(stderr, "cannot open %s\n", filename);
    return -1;
  }

  //  fprintf( stderr, "Reading Gadget file blksize= %d\n", blksize);
  fread( &blksize, 4, 1, fin);
  fread( &gadget_header, sizeof(gadget_header), 1, fin);
  fread( &blksize, 4, 1, fin);
  //  fprintf( stderr, "End Reading Gadget file blksize= %d\n", blksize);

  npart = gadget_header.Npart[1];

  for (i = 0; i < npart; i++) {
    ptcl[i].mass = gadget_header.Massarr[1];
  }

  cache = (float *)malloc(sizeof(float) * 3 * IO_CACHE_SIZE);
  fread( &blksize, 4, 1, fin);
  //  fprintf( stderr, "Reading Gadget file blksize= %d\n", blksize);
  for( int i=0; i<npart; i+=nmemory){
    int nread = nmemory;
    if( (i+nmemory) > npart)  nread = npart - i;
    fread( cache, sizeof(float), 3*nread, fin);
#pragma omp parallel for schedule(auto)
    for( int j=0; j<nread; j++){
      ptcl[i+j].r[0] = cache[3*j+0];
      ptcl[i+j].r[1] = cache[3*j+1];
      ptcl[i+j].r[2] = cache[3*j+2];
    }
  }
  fread( &blksize, 4, 1, fin);
  //  fprintf( stderr, "End Reading Gadget file blksize= %d\n", blksize);

  free(cache);
  fclose(fin);

  return npart;
}

void calc_mesh_density(double *mesh_density, Particle *ptcl,
                       int npart, int ngrid)
{

  for (int p = 0; p < npart; p++) {
    double xt1;
    double dx1;
    double wi11, wi21, wi31, wi12, wi22, wi32, wi13, wi23, wi33;
    int iw11, iw21, iw31, iw12, iw22, iw32, iw13, iw23, iw33;

    xt1  = (double)ptcl[p].r[0] / BOXSIZE * (double)ngrid;
    iw21 = (int)(xt1 + 0.5);
    dx1  = xt1 - (double)iw21;
    wi11 = 0.5 * (0.5 - dx1) * (0.5 - dx1);
    wi21 = 0.75 - dx1 * dx1;
    wi31 = 0.5 * (0.5 + dx1) * (0.5 + dx1);
    iw11 = iw21 - 1;
    iw31 = iw21 + 1;
    if (iw21 == 0) {
      iw11 = ngrid - 1;
    } else if (iw21 == ngrid - 1) {
      iw31 = 0;
    } else if (iw21 == ngrid) {
      iw21 = 0;
      iw31 = 1;
    }

    xt1  = (double)ptcl[p].r[1] / BOXSIZE * (double)ngrid;
    iw22 = (int)(xt1 + 0.5);
    dx1  = xt1 - (double)iw22;
    wi12 = 0.5 * (0.5 - dx1) * (0.5 - dx1);
    wi22 = 0.75 - dx1 * dx1;
    wi32 = 0.5 * (0.5 + dx1) * (0.5 + dx1);
    iw12 = iw22 - 1;
    iw32 = iw22 + 1;
    if (iw22 == 0) {
      iw12 = ngrid - 1;
    } else if (iw22 == ngrid - 1) {
      iw32 = 0;
    } else if (iw22 == ngrid) {
      iw22 = 0;
      iw32 = 1;
    }

    xt1  = (double)ptcl[p].r[2] / BOXSIZE * (double)ngrid;
    iw23 = (int)(xt1 + 0.5);
    dx1  = xt1 - (double)iw23;
    wi13 = 0.5 * (0.5 - dx1) * (0.5 - dx1);
    wi23 = 0.75 - dx1 * dx1;
    wi33 = 0.5 * (0.5 + dx1) * (0.5 + dx1);
    iw13 = iw23 - 1;
    iw33 = iw23 + 1;
    if (iw23 == 0) {
      iw13 = ngrid - 1;
    } else if (iw23 == ngrid - 1) {
      iw33 = 0;
    } else if (iw23 == ngrid) {
      iw23 = 0;
      iw33 = 1;
    }

    wi11 *= ptcl[p].mass;
    wi21 *= ptcl[p].mass;
    wi31 *= ptcl[p].mass;

    RHO(iw11,iw12,iw13)=RHO(iw11,iw12,iw13)+wi11*wi12*wi13;
    RHO(iw21,iw12,iw13)=RHO(iw21,iw12,iw13)+wi21*wi12*wi13;
    RHO(iw31,iw12,iw13)=RHO(iw31,iw12,iw13)+wi31*wi12*wi13;
    RHO(iw11,iw22,iw13)=RHO(iw11,iw22,iw13)+wi11*wi22*wi13;
    RHO(iw21,iw22,iw13)=RHO(iw21,iw22,iw13)+wi21*wi22*wi13;
    RHO(iw31,iw22,iw13)=RHO(iw31,iw22,iw13)+wi31*wi22*wi13;
    RHO(iw11,iw32,iw13)=RHO(iw11,iw32,iw13)+wi11*wi32*wi13;
    RHO(iw21,iw32,iw13)=RHO(iw21,iw32,iw13)+wi21*wi32*wi13;
    RHO(iw31,iw32,iw13)=RHO(iw31,iw32,iw13)+wi31*wi32*wi13;

    RHO(iw11,iw12,iw23)=RHO(iw11,iw12,iw23)+wi11*wi12*wi23;
    RHO(iw21,iw12,iw23)=RHO(iw21,iw12,iw23)+wi21*wi12*wi23;
    RHO(iw31,iw12,iw23)=RHO(iw31,iw12,iw23)+wi31*wi12*wi23;
    RHO(iw11,iw22,iw23)=RHO(iw11,iw22,iw23)+wi11*wi22*wi23;
    RHO(iw21,iw22,iw23)=RHO(iw21,iw22,iw23)+wi21*wi22*wi23;
    RHO(iw31,iw22,iw23)=RHO(iw31,iw22,iw23)+wi31*wi22*wi23;
    RHO(iw11,iw32,iw23)=RHO(iw11,iw32,iw23)+wi11*wi32*wi23;
    RHO(iw21,iw32,iw23)=RHO(iw21,iw32,iw23)+wi21*wi32*wi23;
    RHO(iw31,iw32,iw23)=RHO(iw31,iw32,iw23)+wi31*wi32*wi23;

    RHO(iw11,iw12,iw33)=RHO(iw11,iw12,iw33)+wi11*wi12*wi33;
    RHO(iw21,iw12,iw33)=RHO(iw21,iw12,iw33)+wi21*wi12*wi33;
    RHO(iw31,iw12,iw33)=RHO(iw31,iw12,iw33)+wi31*wi12*wi33;
    RHO(iw11,iw22,iw33)=RHO(iw11,iw22,iw33)+wi11*wi22*wi33;
    RHO(iw21,iw22,iw33)=RHO(iw21,iw22,iw33)+wi21*wi22*wi33;
    RHO(iw31,iw22,iw33)=RHO(iw31,iw22,iw33)+wi31*wi22*wi33;
    RHO(iw11,iw32,iw33)=RHO(iw11,iw32,iw33)+wi11*wi32*wi33;
    RHO(iw21,iw32,iw33)=RHO(iw21,iw32,iw33)+wi21*wi32*wi33;
    RHO(iw31,iw32,iw33)=RHO(iw31,iw32,iw33)+wi31*wi32*wi33;
  }
}

void calc_delta(double *delta, int ngrid)
{
  int64_t datasize = ngrid * ngrid * ngrid;
  double rho_bar;

  rho_bar = 0.0;
#pragma omp parallel for schedule(auto) collapse(3) reduction(+:rho_bar)
  for(int ix=0;ix<ngrid;ix++) {
    for(int iy=0;iy<ngrid;iy++) {
      for(int iz=0;iz<ngrid;iz++) {
        int64_t im = iz + (ngrid+2)*(iy + ngrid*ix);
        rho_bar += delta[im];
      }
    }
  }
  rho_bar /= (double)(datasize);

#pragma omp parallel for schedule(auto) collapse(3)
  for(int ix=0;ix<ngrid;ix++) {
    for(int iy=0;iy<ngrid;iy++) {
      for(int iz=0;iz<ngrid;iz++) {
        int64_t im = iz + (ngrid+2)*(iy + ngrid*ix);
        delta[im] = (delta[im] - rho_bar) / rho_bar;
      }
    }
  }
}

#define cmplx_re(c) ((c)[0])
#define cmplx_im(c) ((c)[1])
#define TINY (1.0e-10)

float correction(float k, float k_N)
{
  float ss = sin(0.5 * PI * k / k_N);
  return (1.0 - SQR(ss) + 0.13333333 * QUART(ss));
}

#define DELTAHAT(a,b,c) (delta_hat[(c) + (ngrid/2+1) * ((b) + ngrid * (a))])

void calc_power(double *delta, int ngrid,
                double **pk, double **kwave, int *nk)
{
  fftw_complex *delta_hat;
  fftw_plan plan;

  float dk;
  double *power, *weight;
  int nkbin;

  plan = fftw_plan_dft_r2c_3d(ngrid, ngrid, ngrid,
                              delta, (fftw_complex *)delta, FFTW_ESTIMATE);

  fftw_execute(plan);

  delta_hat = (fftw_complex *)delta;

  nkbin = ngrid;

  power = (double *)malloc(sizeof(double) * nkbin);
  weight = (double *)malloc(sizeof(double) * nkbin);

  dk = 2.0 * PI / BOXSIZE;

  for (int ik = 0; ik < nkbin; ik++) {
    power[ik] = 0.0;
    weight[ik] = 0.0;
  }

  for (int ix = 0; ix < ngrid; ix++) {
    float xk;
    if (ix <= ngrid / 2) {
      xk = (float)ix;
    } else {
      xk = (float)(ngrid - ix);
    }

    for (int iy = 0; iy < ngrid; iy++) {
      float yk;
      if (iy <= ngrid / 2) {
        yk = (float)iy;
      } else {
        yk = (float)(ngrid - iy);
      }

      for (int iz = 0; iz < ngrid / 2 + 1; iz++) {
        float zk;
        int ik;

        zk = (float)iz;
        ik = (int)(sqrt(SQR(xk) + SQR(yk) + SQR(zk)));

        if (ik < nkbin) {
          float re = cmplx_re(DELTAHAT(ix, iy, iz));
          float im = cmplx_im(DELTAHAT(ix, iy, iz));
          power[ik] += SQR(re) + SQR(im);
          weight[ik] += 1.0;
        }
      }
    }
  }

  {
    float nmesh6 = SQR(CUBE((float)ngrid));
    float k_N = (float)(nkbin / 2) * dk;

    for (int ik = 0; ik < nkbin; ik++) {
      power[ik] /= (weight[ik] + TINY);
      power[ik] /= correction(ik * dk, k_N);
      power[ik] /= nmesh6;
    }
  }

  *nk = nkbin;
  *pk = (double *)malloc(sizeof(double) * nkbin);
  *kwave = (double *)malloc(sizeof(double) * nkbin);

  for (int ik = 0; ik < nkbin; ik++) (*pk)[ik] = power[ik];
  for (int ik = 0; ik < nkbin; ik++) (*kwave)[ik] = dk * ((float)ik + 0.5);

  free(power);
  free(weight);
  fftw_destroy_plan(plan);
}

void calc_potential(double *delta_potential, int ngrid)
{
  fftw_complex *delta_potential_hat;
  fftw_plan forward_plan, backward_plan;
  double dx = BOXSIZE / (double)ngrid;

  forward_plan = fftw_plan_dft_r2c_3d(ngrid, ngrid, ngrid,
                                      delta_potential, (fftw_complex *)delta_potential,
                                      FFTW_ESTIMATE);
  fftw_execute(forward_plan);

  delta_potential_hat = (fftw_complex *)delta_potential;

  for (int ix = 0; ix < ngrid; ix++) {
    int kx;
    if (ix <= ngrid / 2) {
      kx = ix;
    } else {
      kx = ix - ngrid;
    }

    for (int iy = 0; iy < ngrid; iy++) {
      int ky;
      if (iy <= ngrid / 2) {
        ky = iy;
      } else {
        ky = iy - ngrid;
      }

      for (int iz = 0; iz < ngrid / 2 + 1; iz++) {
        int kz = iz;
        int im = iz + (ngrid / 2 + 1) * (iy + ngrid * ix);

        double sx = sin(PI * (double)kx / (double)ngrid);
        double sy = sin(PI * (double)ky / (double)ngrid);
        double sz = sin(PI * (double)kz / (double)ngrid);
        double denom = SQR(sx) + SQR(sy) + SQR(sz);   //分母

        if (denom == 0.0) {
          delta_potential_hat[im][0] = 0.0;
          delta_potential_hat[im][1] = 0.0;
        } else {
          double green = -PI * GRAVITY_G * dx * dx / denom;
          delta_potential_hat[im][0] *= green;
          delta_potential_hat[im][1] *= green;
        }
      }
    }
  }

  backward_plan = fftw_plan_dft_c2r_3d(ngrid, ngrid, ngrid,
                                       delta_potential_hat, delta_potential,
                                       FFTW_ESTIMATE);
  fftw_execute(backward_plan);

#pragma omp parallel for schedule(auto) collapse(3)
  for(int ix = 0;ix < ngrid;ix++) {
    for(int iy = 0;iy < ngrid;iy++) {
      for(int iz = 0;iz < ngrid;iz++) {
        int64_t im = iz + (ngrid + 2) * (iy + ngrid * ix);
        delta_potential[im] /= (ngrid * ngrid * ngrid);
      }
    }
  }

  fftw_destroy_plan(forward_plan);
  fftw_destroy_plan(backward_plan);
}

int output_potential_binary(const char *filename, double *delta_potential, int ngrid)
{
  FILE *fp;

  fp = fopen(filename, "wb");
  if (fp == NULL) {
    fprintf(stderr, "cannot open %s\n", filename);
    return -1;
  }

  for(int ix = 0;ix < ngrid;ix++) {
    for(int iy = 0;iy < ngrid;iy++) {
      int64_t im = (ngrid + 2)*(iy + ngrid * ix);
      if (fwrite(&delta_potential[im], sizeof(double), ngrid, fp) != (size_t)ngrid) {
        fprintf(stderr, "failed to write %s\n", filename);
        fclose(fp);
        return -1;
      }
    }
  }

  fclose(fp);
  return 0;
}


#define __SNAPSHOT_PREFIX__ "../../snapdir_%03d/U2000_%03d_samp0p005.gad.%d"

int output_potential_slice(const char *filename, double *delta_potential,
                           int ngrid, int iz_slice)
{
  FILE *fp;
  double dx = BOXSIZE / (double)ngrid;

  fp = fopen(filename, "w");
  if (fp == NULL) {
    fprintf(stderr, "cannot open %s\n", filename);
    return -1;
  }

  for (int ix = 0; ix < ngrid; ix++) {
    for (int iy = 0; iy < ngrid; iy++) {
      int64_t im = iz_slice + (ngrid + 2) * (iy + ngrid * ix);

      fprintf(fp, "%23.16e %23.16e %23.16e\n",
              dx * ((double)ix + 0.5),
              dx * ((double)iy + 0.5),
              delta_potential[im]);
    }
    fprintf(fp, "\n");
  }

  fclose(fp);
  return 0;
}

int output_delta_slice(const char *filename, double *delta,
                       int ngrid, int iz_slice)
{
  FILE *fp;
  double dx = BOXSIZE / (double)ngrid;

  fp = fopen(filename, "w");
  if (fp == NULL) {
    fprintf(stderr, "cannot open %s\n", filename);
    return -1;
  }

  for (int ix = 0; ix < ngrid; ix++) {
    for (int iy = 0; iy < ngrid; iy++) {
      int64_t im = iz_slice + (ngrid + 2) * (iy + ngrid * ix);

      fprintf(fp, "%23.16e %23.16e %23.16e\n",
              dx * ((double)ix + 0.5),
              dx * ((double)iy + 0.5),
              delta[im]);
    }
    fprintf(fp, "\n");
  }

  fclose(fp);
  return 0;
}



int main(int argc, char **argv)
{
  FILE *fout;
  Particle *ptcl;
  char filename[256];
  char potential_filename[256];
  char delta_slice_filename[256];
  char potential_slice_filename[256];
  char pk_filename[256];
  double *delta;
  double *delta_potential;
  double *pk, *kwave;
  int nkbin;

  if (argc != 2 && argc != 3) {
    fprintf(stderr, "Usage: %s snapshot_index \n", argv[0]);
    return 1;
  }

  int snap_indx = atoi(argv[1]);

  int datasize = NGRID * NGRID * (NGRID+2);
  delta = (double *)malloc(sizeof(double) * datasize);

  for (int i = 0; i < datasize; i++) {
    delta[i] = 0.e0;
  }

  int64_t npart_total = 0;
  for (int32_t i = 0; i < 200; i++) {
    int32_t npart;

    fprintf(stdout, "# Reading %d-th file out of 200.\n", i);

    //    sprintf(filename, "%s.%d", argv[1], i);
    sprintf(filename, __SNAPSHOT_PREFIX__,snap_indx, snap_indx, i);
    printf("%s\n", filename);

    npart = get_gadget_npart(filename);
    if (npart < 0) break;

    printf("%s : npart = %d\n", filename, npart);
    npart_total += npart; printf("%lld\n", npart_total);

    ptcl = (Particle *) malloc(sizeof(Particle)*npart);
    read_gadget_ptcl(filename, ptcl);

    calc_mesh_density(delta, ptcl, npart, NGRID);
    free(ptcl);
  }

  printf("# total number of particles for this snapshot : %lld \n",
         npart_total);

  calc_delta(delta, NGRID);
  sprintf(delta_slice_filename, "delta_slice_%03d.dat", snap_indx);
  output_delta_slice(delta_slice_filename,
                     delta_potential, NGRID, NGRID / 2 );

  delta_potential = (double *)malloc(sizeof(double) * datasize);
  for (int i = 0; i < datasize; i++) {
    delta_potential[i] = delta[i];
  }

  calc_power(delta, NGRID, &pk, &kwave, &nkbin);
  calc_potential(delta_potential, NGRID);

  sprintf(potential_slice_filename, "potential_slice_%03d.dat", snap_indx);
  output_potential_slice(potential_slice_filename,
                         delta_potential, NGRID, NGRID / 2);

#if 0
  fout = fopen("delta_tsc.dat2", "w");
  for (int i = 0; i < NGRID; i++) {
    for (int j = 0; j < NGRID; j++) {
      for (int k = 0; k < NGRID; k++) {
        fprintf(fout, "%d %d %d %14.6e\n",
                i, j, k,
                delta[k + NGRID * (j + NGRID * i)]);
      }
    }
  }
#endif

  sprintf(pk_filename, "power_%03d.dat", snap_indx);
  fout = fopen(pk_filename, "w");
  for (int ik = 1; ik < nkbin / 2; ik++) {
    fprintf(fout, "%14.6e %14.6e\n", kwave[ik], pk[ik]);
  }


  free(pk);
  free(kwave);
  free(delta);
  fclose(fout);


  //printf("# total number of particle in this snapshot : %ld\n", npart_total);
  //printf("# total selected particles = %lld\n", total_selected);
  //printf("wrote zslice_xy.dat\n");

  return 0;
}

#undef DELTAHAT
#undef TINY
#undef cmplx_re
#undef cmplx_im
