# HaussCM: Program to calculate the Haussdorf Chirality Measure (HCM)

This tool was produced to calculate the HCM employing a brute-force approach: a fine grid in rotational space is explored in order to discover the smallest possible Haussdorf distance between specular images of an atomic system. The exploration is parallelized across CPU cores using OpenMP.

## Compilation

Compiling HaussCM only requires a Fortan compiler and OpenMP. For example, using the GNU compiler:

```sh
gfortran -fopenmp -O3 -o haussCM file_parser.f90 hcm_core.f90 hcm_main.f90
```

## Usage

To run HaussCM a valid XYZ file with atomic coordinates should be provided. ´cysteine.xyz´ is provided as an example (its HCM value should be 0.1626):

```sh
./haussCM cysteine.xyz
```

### Options

```sh
-h, --help      Show  help message
-xyz            Print xyz best fit coordinates for the enantiomers.
Files produced by -xyz option:
                Centered reference: 'centered_ref.xyz'.
                Best fit obtained enantiomer: 'enan_fit.xyz'.
```

### OpenMP parallelization

The program will use the default number of threads/cores configured for OpenMP runs. To modify this, change the OMP_NUM_THREADS variable. For example, to run with 4 cores:

```sh
OMP_NUM_THREADS=4 ./haussCM cysteine.xyz
```

