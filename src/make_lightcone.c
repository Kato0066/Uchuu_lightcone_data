/*Copyright (c) 2020 Tomoaki Ishiyama
  MIT License*/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>


const double Gadget_UnitLength_in_Mpc = 1.0;            // 1Mpc
const double Gadget_UnitMass_in_Msun = 1.0e10;          // 1e10 Msun
const double Gadget_UnitVelocity_in_cm_per_s = 1e5;     //  1 km/sec

#define IO_CACHE_SIZE 2097152
#define ZMIN 0.0
#define ZMAX 0.1


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



int main(int argc, char **argv)
{
  FILE *fout;
  Particle *ptcl;
  char filename[256];
  int64_t total_selected = 0;

  if (argc != 2) {
    fprintf(stderr, "Usage: %s base_filename\n", argv[0]);
    return 1;
  }

  fout = fopen("zslice_xy.dat", "w");
  if (fout == NULL) {
    fprintf(stderr, "cannot open zslice_xy.dat\n");
    return 1;
  }

  int64_t npart_total = 0;

  for (int32_t i = 0; i < 200; i++) {
    int32_t npart;

    fprintf(stdout, "# Reading %d-th file out of 200.\n", i);

    sprintf(filename, "%s.%d", argv[1], i);

    npart = get_gadget_npart(filename);

    npart_total += npart;
    if (npart < 0) break;

    ptcl = (Particle *) malloc(sizeof(Particle)*npart);
    if (ptcl == NULL) {
      fprintf(stderr, "malloc failed for %s\n", filename);
      fclose(fout);
      return 1;
    }

    read_gadget_ptcl(filename, ptcl);

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

    free(ptcl);
  }

  fclose(fout);

  printf("# total number of particle in this snapshot : %ld\n", npart_total);
  printf("# total selected particles = %lld\n", total_selected);
  printf("wrote zslice_xy.dat\n");

  return 0;
}
