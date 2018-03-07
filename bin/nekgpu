#!/bin/bash

set -e

DBG=
# Correct format for qsub
# REQ_WALLTIME="00:30:00"
# Correct format for bsub
REQ_WALLTIME="10"
PROJECT="CSC262"

if [ $2 -gt 1024 ]; then
PROSIZE="prod"
else
PROSIZE="prod-devel"
fi
SUBMIT_ARGS=""
WAIT="0"
while true; do
  case "$1" in
    -h|--help )
      echo "Usage: $0 [options] [.rea stem]"
      echo
      echo "Usable options:"
      echo "-h --help: Get help"
      echo "-d --debug: Run in debugger"
      echo "-n|--nodes N: Set number of nodes to N"
      echo "-t NN:NN:NN: requested amount of computer time"
      echo "-s|--submit-args \"-arg1 -arg2\" extra arguments to qsub"
      echo "-w|--wait wait until job is completed"
      exit 1
      shift
      ;;
    -d|--debug )
      echo "*** running in debugger"
      DBG="gdb"
      shift
      ;;
    -n|-nodes|--nodes )
      shift
      CORECOUNT="$1"
      shift
      ;;
    -t )
      shift
      REQ_WALLTIME="$1"
      shift
      ;;
    -s|--submit-args )
      shift
      SUBMIT_ARGS="$1"
      shift
      ;;
    -w|--wait )
      shift
      WAIT="1"
      ;;
    * )
      break
      ;;
  esac
done

rm -f xxt_map.rea

# automatically find .rea file, if unique
if test "$1" = ""; then
  COUNTREA=`ls *.rea | wc -l`
  if test $COUNTREA = 1; then
    REAFILE=`ls *.rea`
    echo "*** found only $REAFILE, picking that one"
  else
    echo "Must specify .rea file; there is more than one here:"
    ls *.rea | cat
    exit 1
  fi
else
  REAFILE=$1
  CORECOUNT=$2
  NODECOUNT=$3
fi

SESSION=${REAFILE%.rea}

if test -d /bgsys; then  # running on BG/Q

  if test -d /gpfs; then # TODO: change this when Mira comes online

      if test "$CORECOUNT" = ""; then
	  CORECOUNT=4
	  echo "*** defaulting to $CORECOUNT nodes on bgq"
      fi
      echo "*** running on bgq with $CORECOUNT cores per node "

      rm -f $SESSION.output
      rm -f logfile
      rm -f xxt_map.rea

      OUTFILE="`pwd`/$SESSION.np=$CORECOUNT-bgq-`date "+%F_%H_%M_%S"`"
      touch $SESSION.rea
      touch $OUTFILE.output
      ln $OUTFILE.output $SESSION.output
      ln $OUTFILE.output logfile

      rm -Rf vtk
      mkdir -p vtk

      echo "qsub  -node $3  --mode c$CORECOUNT -A $PROJECT -t $REQ_WALLTIME -O $OUTFILE nekcem $SESSION"
      COBALTJOB=`qsub -n $3  --mode c$CORECOUNT -A $PROJECT -t $REQ_WALLTIME -O $OUTFILE nekcem $SESSION`
      echo "=== cobalt job $COBALTJOB submitted to veas"

  else

      if test "$CORECOUNT" = ""; then
	  CORECOUNT=4
	  echo "*** defaulting to $CORECOUNT nodes on bgp"
      fi
      echo "*** running on bgp with $CORECOUNT nodes"

      rm -f $SESSION.output
      rm -f logfile
      rm -f xxt_map.rea

      OUTFILE="`pwd`/$SESSION.np=$CORECOUNT-bgsys-`date "+%F_%H_%M_%S"`"
      touch $SESSION.rea
      touch $OUTFILE.output
      ln $OUTFILE.output $SESSION.output
      ln $OUTFILE.output logfile

      rm -Rf vtk
      mkdir -p vtk

      echo "cqsub  -n $CORECOUNT -m vn -p $PROJECT -q $PROSIZE -e BG_MAPPING=TXYZ -t $REQ_WALLTIME -O $OUTFILE nekcem $SESSION"
      COBALTJOB=`cqsub -n $CORECOUNT -m vn -p $PROJECT -q $PROSIZE -e BG_MAPPING=TXYZ -t $REQ_WALLTIME -O $OUTFILE nekcem $SESSION`
      echo "=== cobalt job $COBALTJOB submitted"

  fi

  if test "$WAIT" = "1"; then
    echo "... waiting for job, step 1 "
    zinfo -c $COBALTJOB -w > /dev/null || true
    echo "... waiting for job, step 2"
    zinfo -c $COBALTJOB -e > /dev/null || true
    echo "... waiting for job, step 3"
    while cqstat | grep $COBALTJOB > /dev/null; do
      sleep 1
    done
    echo "--------------------"
    echo "last_error contains:"
    echo "--------------------"
    cat last_error
    echo "--------------------"
    echo "last_output contains:"
    echo "--------------------"
    cat last_output
    echo "=== job finished"
  fi

elif [[ $HOSTNAME =~ .*theta.* ]]; then
    QSUBSCRIPT=$(mktemp -p $PWD)
    echo "#!/bin/bash" >> $QSUBSCRIPT
    echo "aprun -n $CORECOUNT -d 2 -cc depth -j 1 ./nekcem $SESSION" >> $QSUBSCRIPT
    chmod a+x $QSUBSCRIPT
    qsub -A $PROJECT -t $REQ_WALLTIME -q default -n $3 -o logfile $QSUBSCRIPT
    echo "job submitted on ALCF Theta, #CPUs=$CORECOUNT, #nodes=$NODECOUNT"

elif [[ $HOSTNAME =~ .*titan.* ]]; then
    QSUBSCRIPT=$(mktemp -p $PWD)
    DATE=$(date "+%F_%H_%M_%S")
    PROCSPERNODE=$(($CORECOUNT/$NODECOUNT))

    echo "#!/bin/bash -l" >> $QSUBSCRIPT
    echo "#PBS -A $PROJECT" >> $QSUBSCRIPT
    echo "#PBS -N $SESSION" >> $QSUBSCRIPT
    echo "#PBS -o $PWD/$SESSION.np=$CORECOUNT-titan-gpu-$DATE.output" >> $QSUBSCRIPT
    echo "#PBS -e $PWD/$SESSION.np=$CORECOUNT-titan-gpu-$DATE.error" >> $QSUBSCRIPT
    echo "#PBS -l walltime=$REQ_WALLTIME,nodes=$3" >> $QSUBSCRIPT
    echo "#PBS -j oe" >> $QSUBSCRIPT
    echo " cd `pwd`">> $QSUBSCRIPT
    echo " export MPICH_RDMA_ENABLED_CUDA=1" >> $QSUBSCRIPT
    echo " export LD_LIBRARY_PATH=$CRAY_LD_LIBRARY_PATH:$LD_LIBRARY_PATH" >> $QSUBSCRIPT
    echo " CUDA_PROFILE=1 ">> $QSUBSCRIPT
    echo " aprun -n $CORECOUNT -N 1 ./nekcem $SESSION">> $QSUBSCRIPT

    qsub -V $QSUBSCRIPT
    rm $QSUBSCRIPT
    echo "job submitted on OLCF Titan, #CPUs=$CORECOUNT, #Nodes=$NODECOUNT"

elif [[ $HOSTNAME =~ .*summitdev.* ]]; then
    BSUBSCRIPT=$(mktemp -p $PWD)
    DATE=$(date "+%F_%H_%M_%S")
    PROCSPERNODE=$(($CORECOUNT/$NODECOUNT))
 
    #uncomment the following line if GPUDIRECT=true
    #source $OLCF_SPECTRUM_MPI_ROOT/jsm_pmix/bin/export_smpi_env -gpu

    echo "#!/bin/bash" >> $BSUBSCRIPT
    echo "#BSUB -P $PROJECT" >> $BSUBSCRIPT
    echo "#BSUB -J $SESSION" >> $BSUBSCRIPT
    echo "#BSUB -o $PWD/$SESSION.np=$CORECOUNT-summitdev-gpu-$DATE.output" >> $BSUBSCRIPT
    echo "#BSUB -e $PWD/$SESSION.np=$CORECOUNT-summitdev-gpu-$DATE.error" >> $BSUBSCRIPT
    echo "#BSUB -W $REQ_WALLTIME" >> $BSUBSCRIPT
    echo "#BSUB -nnodes $NODECOUNT" >> $BSUBSCRIPT
    echo "export PGI_ACC_NOTIFY=1" >> $BSUBSCRIPT 
    echo "jsrun -n$CORECOUNT -a1 -g1 ./nekcem $SESSION" >> $BSUBSCRIPT

    bsub $BSUBSCRIPT       
    rm $BSUBSCRIPT
    echo "job submitted on OLCF Summitdev, #GPUs=$CORECOUNT, #nodes=$NODECOUNT"

else
    if test "$2" = ""; then
	echo "This is to run with MPI: must specify np# $2"
	exit 1
    fi

    echo "Job to be submitted with np=$2 $SESSION"
    mpiexec -np $2 ./nekcem $SESSION > $SESSION.np=$2.output
fi



