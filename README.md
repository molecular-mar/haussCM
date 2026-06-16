# HaussCM: Program to calculate the Haussdorf Chirality Measure (HCM)

## Compilation

Compiling HaussCM only requires a Fortan compiler and OpenMP. For example, using the GNU compiler:

```sh
gfortran -fopenmp -O3 -o haussCM file_parser.f90 hcm_core.f90 hcm_main.f90
```

## Usage

To run HaussCM a valid XYZ file with atomic coordinates should be provided. ´cysteine.xyz´ is provided as an example:

```sh
./haussCM cysteine.xyz
```

The program will use the default number of threads/cores configured for OpenMP runs. To modify this, change the OMP_NUM_THREADS variable. For example, to run with 4 cores:

```sh
OMP_NUM_THREADS=4 ./haussCM cysteine.xyz
```

