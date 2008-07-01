#! /bin/sh
#
# Automatically install GNU Fortran 4.1.x with dependencies.
#
# Andreas Kloeckner <andreas@mcs.anl.gov>, 8/11/2006
# This script is capable of resuming its work, so you may run it in
# the same directory without deleting intermediate results.

set -e

if test "$1" = ""; then
  echo "usage: $0 INSTALLDIR"
  exit 1
fi
PFX=$1
GCCSUFFIX=-4.1

GMPTRUNK=gmp-4.2.1
MPFRTRUNK=mpfr-2.2.0
GCCTRUNK=gcc-4.1.1
GMPBALL=$GMPTRUNK.tar.bz2
MPFRBALL=$MPFRTRUNK.tar.bz2
GCCBALL=$GCCTRUNK.tar.bz2

mkdir -p dl
cd dl
if ! test -f $GMPBALL; then
  wget -c ftp://ftp.gnu.org/gnu/gmp/$GMPBALL
fi
if ! test -f $MPFRBALL; then
  wget -c http://www.mpfr.org/mpfr-current/$MPFRBALL
fi
if ! test -f $GCCBALL; then
  wget -c ftp://ftp.fu-berlin.de/unix/languages/gcc/releases/$GCCTRUNK/$GCCBALL
fi
cd ..

mkdir -p build
cd build
if ! test -f gmp-done; then
  rm -Rf $GMPTRUNK
  tar xvfj ../dl/$GMPBALL
  cd $GMPTRUNK
  ./configure --prefix=$PFX --disable-static
  make
  make check
  make install
  cd ..
  touch gmp-done
fi

if ! test -f mpfr-done; then
  rm -Rf $MPFRTRUNK
  tar xvfj ../dl/$MPFRBALL
  cd $MPFRTRUNK
  ./configure --prefix=$PFX --disable-static --enable-shared \
    --with-gmp-include=$PFX/include \
    --with-gmp-lib=$PFX/lib
  make
  make check
  make install
  cd ..
  touch mpfr-done
fi

if ! test -f gcc-done; then
  rm -Rf $GCCTRUNK
  tar xvfj ../dl/$GCCBALL
  cd $GCCTRUNK
  ./configure --prefix=$PFX --enable-languages=fortran \
    --with-mpfr=$PFX --with-gmp=$PFX --program-suffix=-ak-4.1
  make
  make install
  cd ..
  touch gcc-done
fi

echo "All done, congratulations. :)"
