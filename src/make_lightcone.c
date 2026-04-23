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


static int read_block_header(FILE *fin, int *blksize)
{
  if (fread(blksize, sizeof(int), 1, fin) != 1) {
    return -1;
  }
  return 0;
}

static int read_block_footer(FILE *fin, int expected_blksize)
{
  int blksize;

  if (fread(&blksize, sizeof(int), 1, fin) != 1) {
    return -1;
  }
  if (blksize != expected_blksize) {
    fprintf(stderr, "block size mismatch: %d != %d\n", blksize, expected_blksize);
    return -1;
  }
  return 0;
}

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

  if (read_block_header(fin, &blksize) != 0) {
    fclose(fin);
    return -1;
  }

  if (fread(&gadget_header, sizeof(GadgetHeader), 1, fin) != 1) {
    fclose(fin);
    return -1;
  }

  if (read_block_footer(fin, blksize) != 0) {
    fclose(fin);
    return -1;
  }

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

  fin = fopen(filename, "rb");
  if (fin == NULL) {
    fprintf(stderr, "cannot open %s\n", filename);
    return -1;
  }

  if (read_block_header(fin, &blksize) != 0) {
    fclose(fin);
    return -1;
  }
  if (fread(&gadget_header, sizeof(GadgetHeader), 1, fin) != 1) {
    fclose(fin);
    return -1;
  }
  if (read_block_footer(fin, blksize) != 0) {
    fclose(fin);
    return -1;
  }

  npart = gadget_header.Npart[1];

  cache = (float *)malloc(sizeof(float) * 3 * IO_CACHE_SIZE);
  if (cache == NULL) {
    fprintf(stderr, "malloc failed for cache\n");
    fclose(fin);
    return -1;
  }

  if (read_block_header(fin, &blksize) != 0) {
    free(cache);
    fclose(fin);
    return -1;
  }

  for (i = 0; i < npart; i += IO_CACHE_SIZE) {
    int nread = IO_CACHE_SIZE;
    int j;

    if (i + nread > npart) {
      nread = npart - i;
    }

    if (fread(cache, sizeof(float), 3 * nread, fin) != (size_t)(3 * nread)) {
      fprintf(stderr, "failed to read positions from %s\n", filename);
      free(cache);
      fclose(fin);
      return -1;
    }

    for (j = 0; j < nread; j++) {
      ptcl[i + j].r[0] = cache[3 * j + 0];
      ptcl[i + j].r[1] = cache[3 * j + 1];
      ptcl[i + j].r[2] = cache[3 * j + 2];
    }
  }

  if (read_block_footer(fin, blksize) != 0) {
    free(cache);
    fclose(fin);
    return -1;
  }

  free(cache);
  fclose(fin);
  return npart;
}



int main(int argc, char **argv)
{
  FILE *fout;
  Particle *ptcl;
  char filename[256];
  int32_t i;
  long long total_selected = 0;

  if (argc != 2) {
    fprintf(stderr, "Usage: %s base_filename\n", argv[0]);
    return 1;
  }

  fout = fopen("zslice_xy.dat", "w");
  if (fout == NULL) {
    fprintf(stderr, "cannot open zslice_xy.dat\n");
    return 1;
  }


  for (i = 0; i < 200; i++) {
    int32_t npart;
    int32_t j;
    int32_t nsel = 0;

    sprintf(filename, "%s.%d", argv[1], i);

    npart = get_gadget_npart(filename);
    if (npart < 0) break;

    ptcl = (Particle *) malloc(sizeof(Particle)*npart);
    if (ptcl == NULL) {
      fprintf(stderr, "malloc failed for %s\n", filename);
      fclose(fout);
      return 1;
    }

    if (read_gadget_ptcl(filename, ptcl) < 0) {
      free(ptcl);
      fclose(fout);
      return 1;
    }

    /*printf("%s : npart = %d\n", filename, npart);
      if (npart > 0) {
      printf("  ptcl[0] r = (%g, %g, %g)\n",
      ptcl[0].r[0],
      ptcl[0].r[1],
      ptcl[0].r[2]);
      }
    */

    for (j = 0; j < npart; j++) {
      if (ptcl[j].r[2] >= ZMIN && ptcl[j].r[2] <= ZMAX) {
        fprintf(fout, "%.10g %.10g %.10g\n",
                ptcl[j].r[0],
                ptcl[j].r[1],
                ptcl[j].r[2]);
        nsel++;
      }
    }

    total_selected += nsel;
    printf("  selected in %.1f <= z <= %.1f : %d\n", ZMIN, ZMAX, nsel);

    free(ptcl);
  }

  fclose(fout);

  printf("total selected particles = %lld\n", total_selected);
  printf("wrote zslice_xy.dat\n");

  return 0;
}
