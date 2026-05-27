CC=gcc
CFLAGS = -O3 -funroll-loops -ffast-math
CFLAGS += -fopenmp

CLIBS=-lm

all: make_lightcone

MAKE_LIGHTCONE_OBJ = make_lightcone.o
MAKE_LIGHTCONE_DEP := $(MAKE_LIGHTCONE_OBJ)

FFTW_ROOT = /opt/fftw-3.3.10
CFLAGS += -I$(FFTW_ROOT)/include
FFTW_LIBS = -L$(FFTW_ROOT)/lib -lfftw3


make_lightcone: $(MAKE_LIGHTCONE_DEP)
	$(CC) $(CFLAGS) -o $@ $(MAKE_LIGHTCONE_OBJ) $(FFTW_LIBS) $(CLIBS)

clean:
	rm -rf *.o

distclean: clean
	rm -rf *~

run:
#	./make_lightcone ../../snapdir_050/U2000_050_samp0p005.gad > result.txt
	./make_lightcone /import/cage1/kato/Uchuu/snapdir_050/U2000_050_samp0p005.gad | tee result.txt
