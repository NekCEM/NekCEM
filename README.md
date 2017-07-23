# NekCEM

NekCEM is a discontinous-Galerkin, spectral-element solver for
Maxwell's equations and the drift-diffusion equations written in
Fortran and C.  It runs efficiently in parallel on a wide variety
of systems, from laptops to the supercomputers at the Argonne
Leadership Computing Facility (ALCF) and the Oak Ridge Leadership
Computing Facility (OLCF). Its core is based on the computational
fluid dynamics code [Nek5000][Nek5000].

## Installing

### Dependencies

To run simulations with NekCEM you will need the following
things.

- An MPI implementation.
- Python 2.7 or higher (including all versions of Python 3).
- BLAS and Lapack.

Some notes on the dependencies:

- To keep things simple, make sure the compiler wrappers
  `mpif77` (or `mpifort`) and `mpicc` are on your path. This
  isn't strictly necessary, but without them you will have to do
  more work when compiling simulations.
- Python is only used in the build process.
- The system version of Python on some ALCF and OLCF systems is
  2.6; use `softenv` or `modules` to switch to a more recent
  version. Run `soft add +python` on a `softenv` system and
  `module load python` on a modules system.
- Again to keep things simple, make sure you can link to BLAS and
  Lapack using `-lblas` and `-llapack`.

### Standard install

To install NekCEM run

```
git clone https://github.com/NekCEM/NekCEM
cd NekCEM
sudo make install
```

The command `make install` does a couple of things.

- It copies `NekCEM/src` and `NekCEM/bin` to `/usr/local`.
- It symlinks some scripts to `/usr/local/bin`.

Note that installing to `/usr/local` is simply the default option; the
install directory can be controlled in the standard way using the
variables `DESTDIR`, `prefix`, and `bindir`.

### Development install

If you want to help develop the NekCEM source code, first fork the
NekCEM repo on Github. Then do

```
git clone https://github.com/<github-username>/NekCEM
cd NekCEM
git remote add upstream https://github.com/NekCEM/NekCEM
sudo make install_inplace
```

The command `sudo make install_inplace` only symlinks scripts to
`/usr/local/bin`, allowing a developer to edit the source in
their local clone while still having the necessary scripts on
their path.

## Running simulations with NekCEM

Setting up a simulation with NekCEM requires creating four files.

- A user file. This is a Fortran file which contains various
  subroutines used to control the solvers. Its file extension should
  be `usr`.
- A size file. This file contains compile-time parameters. It should
  be called `SIZE`.
- A read file. This file contains parameters which are read at
  runtime. Its file extension should be `rea`.
- A map file. This file contains the mapping between processors and
  elements.  Its file extension should be `map`, and it must have the
  same stem as the read file.

A typical NekCEM simulation will be set up like this

```
example
├── readfile.map
├── readfile.rea
├── userfile.usr
└── SIZE
```

To build and run the code do the following from the `example`
directory.

```
configurenek <solver> userfile
make
mpirun -np <number-of-processors> ./nekcem readfile &> log
```

Let's break down what's going on.

- In the first step `configurenek` creates a makefile. The `<solver>`
  option determines which equations the application is targeting; it
  should be one of `maxwell`, `drift`, or `shrod`.
- In the second step the makefile builds the code in the normal way;
  it produces an executable `nekcem`.
- In the third step the code is run in the normal way for MPI
  applications.

The third step can be replaced with `nek readfile
<number-of-processors>`. On a typical system this will do the exact
same thing as `mpirun`, but on ALCF and OLCF machines it will also
queue your job correctly.

### Running the Tests

The tests can be run with `bin/runtests [options]`. For a complete
list of options use the `-h` flag.

[Nek5000]: https://github.com/Nek5000/Nek5000
