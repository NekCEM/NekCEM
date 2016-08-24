#!/bin/bash -x
#COBALT -t 30
#COBALT -q cache-quad
#COBALT -A GasJetsCyl_tesp

rpn=$1
thr=1
case=box        

echo $case     >  SESSION.NAME
echo `pwd`'/' >>  SESSION.NAME
touch $case.rea
rm -f ioinfo
mv -f $case.his $case.his1
mv -f $case.sch $case.sch1

aprun -n $((COBALT_JOBSIZE*rpn)) \
      -N $rpn \
      -d $thr \
      -cc depth \
      -j 1 \
      ./nekcem   

exit $?


