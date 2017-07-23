# NekCEM

## Installing

### Standard install

To install NekCEM run

```
git clone https://github.com/NekCEM/NekCEM
cd NekCEM
sudo make install
```

The command `make install` does a couple of things.

- It copies `NekCEM/src` and `NekCEM/bin` to `/usr/local`.
- It symlinks the scripts `configurenek` and `nek` to
  `/usr/local/bin`.

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

The command `sudo make install_inplace` only symlinks the scripts
`configurenek` and `nek` to `/usr/local/bin`, allowing a developer to
edit the source in their local clone while still having `configurenek`
and `nek` on their PATH.

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
same thing as `mpirun`, but on a queued system like Mira it will also
queue your job correctly.

### Running the Tests

The tests can be run with `bin/runtests [options]`. For a complete
list of options use the `-h` flag.
