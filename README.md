# NekCEM

## Building and running

Suppose that NekCEM has been installed into the directory `NEK` and
you have the following example.

```
example
|- example.usr
|- example.rea
|- example.box
\- SIZE
```

### Using `makenek`

To build this example:

- copy the script `$NEK/bin/makenek` to `example`
- set the script variable `NEK` to the location of your NekCEM
  installation
- set `APP` to the application you are targeting; should be one of
  `maxwell`, `drift`, or `schrod`
- set `FC` and `CC` to your desired Fortran and C compilers
- run `./makenek example`

This should create an executable `nekcem` in `example`. By default
`makenek` will try to determine a sensible set of `FFLAGS` and
`CFLAGS`, but you can override it by setting `FFLAGS` and `CFLAGS`.
Similarly it will try to determine reasonable values for `LD` and
`LDFLAGS`. You can also add additional flags to the set of flags it
chooses by setting `EXTRAFFLAGS`, `EXTRACFLAGS`, and `EXTRALDFLAGS`.

To run the example run `$NEK/bin/nek example $NP`, where `NP` is the
number of processors you want to use for the run.

### Using `configurenek`

The script `makenek` is a front-end for the script `configurenek`,
which can also be called directly. To build with `configurenek`:

- run `$NEK/bin/configurenek app example [options]` in the `example`
  directory. As with `makenek`, the argument `app` is the application
  being targeted. There are many possible options; run with the `-h`
  flag for a full list. The most basic options are `--FC` and `--CC`,
  which you can use to set the Fortran and C compilers being
  used. Alternatively, if you are on a known machine such as `cetus`,
  then you can just pass `--arch cetus`, which will correctly set the
  compilers and flags for you. See `$NEK/bin/arch.json` for a list of
  known machines.
- Running `configurenek` creates a makefile in `example` which you can
  use to compile the example.

### Running the Tests

To run the tests, run `$NEK/bin/runtests [options]`. For a complete
list of the options use the `-h` flag.
