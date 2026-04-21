CC=gcc
CFLAGS=-O3 -funroll-loops -ffast-math
CFLAGS+= -fopenmp

CLIBS=-lm

all: make_lightcone

MAKE_LIGHTCONE_OBJ = make_lightcone.o
MAKE_LIGHTCONE_DEP := $(MAKE_LIGHTCONE_OBJ)


make_lightcone: $(MAKE_LIGHTCONE_DEP)
	$(CC) $(CFLAGS) -o $@ $(MAKE_LIGHTCONE_OBJ) $(CLIBS)

clean:
	rm -rf *.o

distclean: clean
	rm -rf *~
