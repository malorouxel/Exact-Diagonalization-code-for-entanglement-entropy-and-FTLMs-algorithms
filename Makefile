# Serial build.   make
# MPI build.      make mpi
# Intel MKL.      make MKL=1   or   make mpi MKL=1
# Tests.          make test

CC      = gcc
MPICC   = mpicc
CFLAGS  = -O3 -std=c11 -march=native -fopenmp -Wall -D_GNU_SOURCE
LIBS    = -llapacke -lopenblas -lm

ifeq ($(MKL),1)
CC      = icx
MPICC   = mpiicx
CFLAGS  = -O3 -std=c11 -march=native -qopenmp -qmkl -Wall -D_GNU_SOURCE -DUSE_MKL
LIBS    = -lm
endif

SRC     = parallel.c tables.c indexing.c sector.c hamiltonian.c lanczos.c observables.c entropy.c thermodynamics.c io.c logger.c
SERIAL  = $(addprefix build/serial/, $(SRC:.c=.o))
PARALLEL = $(addprefix build/mpi/, $(SRC:.c=.o))
HEADERS = $(wildcard src/*.h)

all: bin/lanczos bin/entropy bin/thermo

mpi: bin/lanczos_mpi bin/entropy_mpi

bin/lanczos: $(SERIAL) build/serial/main.o | bin
	$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

bin/entropy: $(SERIAL) build/serial/entropy_main.o | bin
	$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

bin/thermo: $(SERIAL) build/serial/thermo_main.o | bin
	$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

bin/lanczos_mpi: $(PARALLEL) build/mpi/main.o | bin
	$(MPICC) $(CFLAGS) -DUSE_MPI $^ -o $@ $(LIBS)

bin/entropy_mpi: $(PARALLEL) build/mpi/entropy_main.o | bin
	$(MPICC) $(CFLAGS) -DUSE_MPI $^ -o $@ $(LIBS)

build/serial/%.o: src/%.c $(HEADERS) | build/serial
	$(CC) $(CFLAGS) -c $< -o $@

build/mpi/%.o: src/%.c $(HEADERS) | build/mpi
	$(MPICC) $(CFLAGS) -DUSE_MPI -c $< -o $@

bin build/serial build/mpi:
	mkdir -p $@

test: all mpi
	bash tests/run_tests.sh

clean:
	rm -rf build bin

.PHONY: all mpi test clean
