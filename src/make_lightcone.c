/*Copyright (c) 2020 Tomoaki Ishiyama
  MIT License*/

#include <stdio.h>
#include <stdlib.h>

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

    fread(&blksize, 4, 1, fin);
    fread(&gadget_header, sizeof(GadgetHeader), 1, fin);
    fread(&blksize, 4, 1, fin);

    fclose(fin);

    return gadget_header.Npart[1];
}

int main(int argc, char **argv)
{
    int npart;

    npart = get_gadget_npart(argv[1]);
    printf("npart = %d\n", npart);

    return 0;
}
