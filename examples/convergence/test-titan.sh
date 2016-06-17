#!/bin/bash

# this is for MPI runs on Vesta/Cetus/Mira

 nproc=32
 sizeFile="SIZEu"

for proc in 'MPI'; do
#   for ele in 4 5 6 7 8 9 10 11 12; do
    for ele in 12  ; do
        eleg=$ele*$ele*$ele

#       for N in 2 4 6 8 10 12 14 16 18; do
        for N in 18 16 14 12 10 8 6 4 2  ; do
	    dir=${proc}_E${ele}_N${N}
             
	    #Make new directory
	    mkdir ../${dir}
	    cd ../${dir}
	    echo $proc

	    #copy in the right files
	    cp ../convergence/SIZEu .
	    cp ../convergence/base.usr 3dbox_E=${ele}_b.usr
	    cp ../convergence/3dbox_E=${ele}_b.* .
	    reaFile="3dbox_E=${ele}_b.rea"

	    #Change rea file
            sed -i '/NSTEPS/c\  1000                         11: x NSTEPS: Total # of timesteps' $reaFile
            sed -i '/Timestep Size/c\  -0.005                       12: x DT:  Timestep Size with (-), eg. -0.05; CF' $reaFile
            sed -i '/IOCOMM/c\  1000                          13: x IOCOMM : Print statement at every IOCOMM' $reaFile
            sed -i '/IOSTEP/c\  1000                          15: x IOSTEP : Produce outputs at every IOSTEP' $reaFile

	    #Change SIZEu file
            sed -i '/polynomial order/c\      parameter (lxi =   '$N')  ! polynomial order' $sizeFile
            sed -i '/number of mpi/c\      parameter (lp  =   '$nproc')  ! number of mpi' $sizeFile 
            sed -i '/number of lelg/c\      parameter (lelg=   '$eleg')  ! number of lelg' $sizeFile

	    #Submit job 
                ../../bin/cleanall   
                ../../bin/makenekmpi -a titan-pgi-mpi
                ../../bin/nek 3dbox_E=${ele}_b 32 2  
                echo "Hello",3dbox_E=${ele}_b,$N,$ele
	done
  done
done
