#!/bin/bash

# this is for MPI runs on Vesta/Cetus/Mira

#nproc=32
 nproc=4 
 sizeFile="SIZEu"

for proc in 'MPI'; do
#   for ele in 3 4 5 6 7 8 9 10 11 12; do
    for ele in 3 4 ; do
        eleg=ele^3

#       for N in 2 4 6 8 10 12 14 16; do
        for N in 2 4 6 ; do
	    dir=${proc}_E${ele}_N${N}
             
	    #Make new directory
	    mkdir ../${dir}
	    cd ../${dir}
	    echo $proc

	    #copy in the right files
	    cp ../convergence/SIZEu .
	    cp ../convergence/base.usr .
	    cp ../convergence/3dbox_E=${ele}_b.* .
	    
	    #Change SIZEu file
            sed -i '/polynomial order/c\      parameter (lxi =   '$N')  ! polynomial order' $sizeFile
            sed -i '/number of mpi/c\      parameter (lp  =   '$nproc')  ! number of mpi' $sizeFile 
            sed -i '/number of lelg/c\      parameter (lelg=   '$ele')  ! number of lelg' $sizeFile

	    #Submit job 
                ../../bin/cleanall   
                ../../bin/makenekmpi -a linux-gnu-mpi
                ../../bin/nek box $nproc
                echo "Hello",$nproc,$N,$ele
	done
  done
done
