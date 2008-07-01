#! /bin/bash

set -e
findbase()
{
  for basedir in . .. ../.. ../../.. ../../../..; do
    if test -d $basedir/src; then
      echo $basedir
      return
    fi
  done
  echo Could not find base directory. Sorry.  > /dev/stderr
  exit 1
}
BASEDIR=`findbase`
echo "base directory is $BASEDIR"

cd $BASEDIR

for i in examples/*; do
  if ls $i/*.usr &> /dev/null ; then
    echo "generating tags in $i"
    (cd $i; etags *.usr SIZE*)
  fi
done

echo "generating tags in src"
(cd src; CAPS_FILES=[A-Z]*; etags *.{F,c} `echo $CAPS_FILES | grep -v TAGS`)
