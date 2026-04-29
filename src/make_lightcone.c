/*Copyright (c) 2020 Tomoaki Ishiyama
  MIT License*/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>


const double Gadget_UnitLength_in_Mpc = 1.0;            // 1Mpc
const double Gadget_UnitMass_in_Msun = 1.0e10;          // 1e10 Msun
const double Gadget_UnitVelocity_in_cm_per_s = 1e5;     //  1 km/sec

#define IO_CACHE_SIZE 2097152
#define NGRID 64
#define RHO(a,b,c) mesh_density[(c) + ngrid * ((b) + ngrid * (a))]
#define BOXSIZE 2000.0
//#define ZMIN 0.0
//#define ZMAX 0.1


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
                       int npart, int ngrid, double boxsize)
{
  int p;

  for (int i = 0; i < ngrid; i++) {
    for (int j = 0; j < ngrid; j++) {
      for (int k = 0; k < ngrid; k++) {
        RHO(i, j, k) = 0.e0;
      }
    }
  }

  for (p = 0; p < npart; p++) {
    double xt1;
    double dx1;
    double wi11, wi21, wi31, wi12, wi22, wi32, wi13, wi23, wi33;
    int iw11, iw21, iw31, iw12, iw22, iw32, iw13, iw23, iw33;

    xt1  = (double)ptcl[p].r[0] / boxsize * (double)ngrid;
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

    xt1  = (double)ptcl[p].r[1] / boxsize * (double)ngrid;
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

    xt1  = (double)ptcl[p].r[2] / boxsize * (double)ngrid;
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



int main(int argc, char **argv)
{
  FILE *fout;
  Particle *ptcl;
  char filename[256];
  double *mesh_density;
  //int64_t total_selected = 0;

  if (argc != 2) {
    fprintf(stderr, "Usage: %s base_filename\n", argv[0]);
    return 1;
  }

  /*fout = fopen("zslice_xy.dat", "w");
    if (fout == NULL) {
    fprintf(stderr, "cannot open zslice_xy.dat\n");
    return 1;
    }


    int64_t npart_total = 0;

    for (int32_t i = 0; i < 200; i++) {
    int32_t npart;

    fprintf(stdout, "# Reading %d-th file out of 200.\n", i);

    sprintf(filename, "%s.%d", argv[1], i);
  */

  int32_t npart;
  npart = get_gadget_npart(argv[1]/*filename*/);

  /*
    npart_total += npart;
    if (npart < 0) break;
  */

  ptcl = (Particle *) malloc(sizeof(Particle)*npart);
  if (ptcl == NULL) {
    fprintf(stderr, "malloc failed for %s\n", filename);
    return 1;
  }

  read_gadget_ptcl(argv[1]/*filename*/, ptcl);

  /*
    int32_t nsel = 0;
    #pragma omp parallel for schedule(auto) reduction(+:nsel)
    for (int32_t j = 0; j < npart; j++) {
    if (ptcl[j].r[2] >= ZMIN && ptcl[j].r[2] <= ZMAX) {
    fprintf(fout, "%14.6e %14.6e %14.6e\n",
    ptcl[j].r[0],
    ptcl[j].r[1],
    ptcl[j].r[2]);
    nsel++;
    }
    }

    total_selected += nsel;
    printf("#  selected in %12.4e <= z <= %12.4e : %d\n", ZMIN, ZMAX, nsel);
  */

  int datasize = NGRID * NGRID * NGRID;
  mesh_density = (double *)malloc(sizeof(double) * datasize);
  if (mesh_density == NULL) {
    fprintf(stderr, "malloc failed for mesh_density\n");
    fclose(fout);
    return 1;
  }

  calc_mesh_density(mesh_density, ptcl, npart, NGRID, BOXSIZE);

  fout = fopen("density_tsc.dat", "w");
  if (fout == NULL) {
    fprintf(stderr, "cannot open density_tsc.dat\n");
    free(mesh_density);
    free(ptcl);
    return 1;
  }

  for (int i = 0; i < NGRID; i++) {
    for (int j = 0; j < NGRID; j++) {
      for (int k = 0; k < NGRID; k++) {
        fprintf(fout, "%d %d %d %14.6e\n",
                i, j, k,
                mesh_density[k + NGRID * (j + NGRID * i)]);
      }
    }
  }

  free(ptcl);
  free(mesh_density);
  fclose(fout);

  printf("npart = %d\n", npart);
  printf("wrote density_tsc.dat\n");

  //printf("# total number of particle in this snapshot : %ld\n", npart_total);
  //printf("# total selected particles = %lld\n", total_selected);
  //printf("wrote zslice_xy.dat\n");

  return 0;
}
