# lanczos-ed

Exact diagonalization and finite temperature Lanczos methods for quantum spin chains.
The code works for spin 1/2 and spin 1 chains with any two-site interactions built from S_i.S_j, (S_i.S_j)^2 and local fields Sz_i.
It runs on one core, with OpenMP threads, or on many nodes with MPI.

The code computes:

- the ground state energy of every magnetization sector
- the ground state vector, saved to disk for later use
- the ground state degeneracy in each sector
- the entanglement entropy of any subsystem, contiguous or not
- the specific heat with three finite temperature Lanczos methods (SFTLM, RFTLM, OFTLM) and jackknife error bars

Small sectors are solved by full exact diagonalization. Large sectors are solved with the Lanczos algorithm.

## Requirements

- a C11 compiler (gcc or Intel icx)
- BLAS and LAPACKE (OpenBLAS or Intel MKL)
- an MPI library for the parallel build (OpenMPI, MPICH or Intel MPI)

On Ubuntu or Debian:

```bash
sudo apt install build-essential libopenblas-dev liblapacke-dev libopenmpi-dev openmpi-bin
```

## Build

```bash
make            # serial build: bin/lanczos bin/entropy bin/thermo
make mpi        # MPI build: bin/lanczos_mpi bin/entropy_mpi
make MKL=1      # use icx and Intel MKL instead of gcc and OpenBLAS
make mpi MKL=1  # use mpiicx and Intel MKL
make test       # build everything and run the tests
make clean
```

The serial and MPI programs are built from the same source files.
The MPI calls are isolated in `src/parallel.c`.

## Quick start

```bash
make
cd examples/aklt_spin_one
../../bin/lanczos
cat AKLT_N8/Results.txt
```

This solves the AKLT chain with 8 sites and periodic boundaries.
The ground state energy is -16/3.

With MPI:

```bash
OMP_NUM_THREADS=1 mpirun -np 4 ../../bin/lanczos_mpi
```

Both programs read `InputsParameters.txt` and `Interactions.txt` in the current folder.
You can give other files on the command line:

```bash
../../bin/lanczos my_parameters.txt my_interactions.txt
```

## Input files

### Parameters

The parameters file alternates a label line and a value line.
The label line is ignored. The order of the entries is fixed.

| Entry | Value |
|---|---|
| Spin | `One` or `Half` |
| Number of sites N | an integer |
| ED threshold | sectors smaller than this are solved by exact diagonalization |
| Sz symmetry | `Yes` to work sector by sector, `No` to work in the full Hilbert space |
| Number of random vectors R | used by the FTLM methods |
| Lanczos iterations | `Conv` to stop when the lowest Ritz value converges, or a fixed number of iterations |
| Ground state degeneracy | `Unknown` or a number. Only used by RFTLM when `dege` is not computed |
| Results | any of `E0 GS dege Entropy SFTLM RFTLM OFTLM` on one line |
| Bipartition size | number L of sites of A for `Entropy`. A is made of the sites 0 to L-1 |
| Output folder | a name. The code appends `_N<N>` |
| Temperatures | `log10(Tmin) log10(Tmax) NT` |
| Seed | optional. `Random` or an integer. Without it the seed comes from the clock |

See `examples/aklt_spin_one/InputsParameters.txt` for a complete file.

### Interactions

The first line is a header. Each other line has the form `J;OP;SITES`.
Sites start at 0.

| OP | Term | SITES |
|---|---|---|
| `SS` | J S_i.S_j | `i-j` |
| `SS2` | J (S_i.S_j)^2 | `i-j` |
| `Sz` | J Sz_i | `i` |

Example for a bilinear-biquadratic chain with periodic boundaries:

```
J;OP;SITES
0.955336;SS;0-1
0.295520;SS2;0-1
0.955336;SS;1-2
0.295520;SS2;1-2
...
0.955336;SS;7-0
0.295520;SS2;7-0
```

Open boundaries are obtained by removing the bond between the last and the first site.

## Output files

All files are written in the output folder.

| File | Content |
|---|---|
| `Results.txt` | lowest energy of each sector, ground state sector and energy, entropy |
| `dege.txt` | ground state degeneracy of each sector |
| `GroundState.bin` | ground state vector |
| `<Sz>_<r>.txt` | Ritz values and weights of random vector r in sector Sz (SFTLM and RFTLM) |
| `O<Sz>_<r>.txt` | same for the Lanczos runs orthogonal to the ground states (OFTLM) |
| `<Sz>_exact.txt` | full spectrum of a sector solved by exact diagonalization |
| `<method>_R<R>_thermo.txt` | T, Z, E/N, C/N, jackknife C/N and the 1 sigma and 2 sigma bounds |
| `Run.log` | duration and memory estimate of each step |

The sector label Sz is the sum of the local states. A local state is 0 or 1 for spin 1/2 and 0, 1 or 2 for spin 1.
Sector Sz has magnetization Sz - N s. With the Sz symmetry only the sectors with a magnetization of 0 or less are computed.
The other sectors are obtained by spin flip.

## Entanglement entropy of a saved ground state

Add `GS` to the Results line to save the ground state.
Then compute the von Neumann entropy of as many subsystems as needed:

```bash
../../bin/entropy AKLT_N8/GroundState.bin 0-3 0,2,4,6 1,5
```

Each argument is one subsystem A. Sites are separated by commas. A range `a-b` adds every site from a to b.
The program prints one line per subsystem:

```
A = 0-3 : S = 1.385377692757438
A = 0,2,4,6 : S = 3.959512508188129
A = 1,5 : S = 2.193520099180299
```

The entropy uses the natural logarithm.
`entropy_mpi` does the same with MPI.

The file `GroundState.bin` starts with a header of 32 bytes: N, spin, Sz symmetry and sector as four 32-bit integers, then the sector dimension as a 64-bit integer and the energy as a double.
The amplitudes follow as doubles in basis order.

## Thermodynamics from a finished run

`thermo` computes the specific heat again from the spectral files.
It can use only the first R random vectors.

```bash
../../bin/thermo InputsParameters.txt 5
```

## How it works

**Basis.** A basis state is a number in base 2 (spin 1/2) or base 3 (spin 1). Each digit is the local state of one site.
With the Sz symmetry the states of one sector are ranked with the combinatorial number system.
A spin 1 state is split into two spin 1/2 states: one marks the sites in state 2, the other gives the state of the remaining sites.

**Lanczos.** The ground state uses two passes. The first pass builds the tridiagonal matrix and finds the lowest Ritz value.
The second pass runs the recursion again from the same random vector and builds the ground state vector.
Only three vectors are stored at any time.

**Degeneracy.** After the ground state, new Lanczos runs are made orthogonal to the known ground states.
A new ground state is added while the lowest Ritz value equals the ground state energy.

**FTLM.** Each random vector gives Ritz values and weights.
SFTLM is the standard method. RFTLM treats the lowest level exactly with its degeneracy.
OFTLM adds the ground states exactly and samples the orthogonal space.
The error bars come from a jackknife over the R random vectors.

**MPI.** Each process owns a contiguous slice of the basis and the same slice of every vector.
Each process computes its own rows of H|v>. Row i needs the values v_j of every state j connected to i.
Before the first product each process lists the indices it needs from other processes. These are the ghost indices.
At each product one `MPI_Alltoallv` sends the ghost values. No process ever stores a full vector.
Random vectors depend only on the seed and the global index. The serial and MPI programs give the same results for the same seed.

**Entropy.** The ground state lives in one Sz sector. Its Schmidt matrix is block diagonal in the charge of the subsystem A.
Each block is diagonalized separately with an SVD. The total memory is the size of the sector instead of the full Hilbert space.
With MPI the blocks are shared between the processes.

## Performance notes

- In the serial program a hash table maps basis states to indices. With MPI the index is computed from the state to save memory.
- The program uses OpenMP in H|v> when `OMP_NUM_THREADS` is above 1. Set `OMP_NUM_THREADS=1` for pure MPI.
- Exact diagonalization stores the full matrix of the sector. Keep the ED threshold at a few thousand states.
- The local slice of a vector must hold fewer than 2^31 elements on each process.

## Repository layout

```
src/
  main.c            ground state, degeneracy, FTLM and entropy driver
  entropy_main.c    entanglement entropy of a saved ground state
  thermo_main.c     thermodynamics from a finished run
  parallel.c        MPI layer, empty in the serial build
  tables.c          binomial coefficients, powers and sector sizes
  indexing.c        ranking of basis states
  sector.c          basis of one sector and MPI ghost exchange
  hamiltonian.c     matrix elements and H|v>
  lanczos.c         Lanczos runs, ground state and degeneracy
  observables.c     exact diagonalization and FTLM files
  entropy.c         entanglement entropy
  thermodynamics.c  specific heat and jackknife
  io.c              input files and ground state file
  logger.c          timing and memory log
examples/           AKLT spin 1 chain and Heisenberg spin 1/2 chain
tests/              checks of the serial and MPI programs
```

## Tests

```bash
make test
```

The tests check the exact AKLT energy, the Heisenberg energy, the agreement between the serial and MPI programs for 1 to 3 processes, and the entropy program.
Set `MPIRUN` to change the MPI launcher, for example `MPIRUN="srun -n" make test`.

## License

MIT. See `LICENSE`.
