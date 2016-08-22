#!/bin/bash

# this is for MPI runs on Vesta/Cetus/Mira

N=7

for proc in 'MPI' ; do
#  for nproc in 1 2 4 8 16 32 64 128 256 512; do
   for nproc in 1 2 ; do
	for ele in 1 2 4 8 16 32 64 128 256 512 1024 2048 4096; do
	    dir=${proc}_${nproc}_${ele}

        if  [ "$ele" -ge "$nproc" ]; then

	    #Make new directory
	    mkdir ../${dir}
	    cd ../${dir}
	    echo $proc

	    #copy in the right files
	    if [ $proc == 'MPI' ]; then
	    	cp ../mpi_scale/SIZEu .
	    	cp ../mpi_scale/box.* .
	       if    [ $ele -ge 32 ]; then
	    	cp ../mpi_scale/data/mesh/b$ele.* .
               fi
	    else
	    	echo 'Must be MPI '
	    	exit
	    fi
	    
	    #Change box file
	    sed -i '/(number of elements/c\-1 -1 -'$ele' (number of elements in x,y,z)' box.box
            sed -i '/polynomial order/c\      parameter (lxi =   '$N')  ! polynomial order' SIZEu
            sed -i '/number of mpi/c\      parameter (lp  =   '$nproc')  ! number of mpi' SIZEu
            sed -i '/number of lelg/c\      parameter (lelg=   '$ele')  ! number of lelg' SIZEu
            sed -i '/number of lelz/c\      parameter (lelx =1,lely =lelx ,lelz =  '$ele')  ! number of lelz' SIZEu
            sed -i '/NEL,NDIM,NELV/c\   '$ele'    -3   '$ele'           NEL,NDIM,NELV' box.rea

            ../../bin/cleanall   
            ../../bin/makenekmpi -a theta-intel

	    #Submit job 
	    if    [ $ele -lt 32 ]; then
	       if    [ $nproc -le 32 ]; then
                    cp ~/run.sh . 
                    cp ~/rr     . 
                    sed -i '/qsub/c\ qsub -n 1 ./run.sh '$nproc'' rr
		    echo 'hello',$proc,$nproc,$ele
               fi
            else   
	       if    [ $nproc == 32 ]; then
                   cp ~/run.sh . 
                    cp ~/rr     . 
                    sed -i '/case/c\case=b'$ele'' run.sh
                    sed -i '/qsub/c\qsub -n 1 ./run.sh '$nproc'' rr
                    echo 'hello',$proc,$nproc,$ele
	       elif  [ $nproc == 64 ]; then
                    cp ~/run.sh . 
                    cp ~/rr     . 
                    sed -i '/case/c\case=b'$ele'' run.sh
                    sed -i '/qsub/c\qsub -n 1 ./run.sh '$nproc'' rr
                    echo 'hello',$proc,$nproc,$ele

               elif  [ $nproc == 128 ]; then
                    cp ~/run.sh . 
                    cp ~/rr     . 
                    sed -i '/case/c\case=b'$ele'' run.sh
                    sed -i '/qsub/c\ qsub -n 2 ./run.sh 64' rr
                    echo 'hello',$proc,$nproc,$ele
               elif  [ $nproc == 256 ]; then
                    cp ~/run.sh . 
                    cp ~/rr     . 
                    sed -i '/case/c\case=b'$ele'' run.sh
                    sed -i '/qsub/c\ qsub -n 4 ./run.sh 64' rr
                echo 'hello',$proc,$nproc,$ele
               elif  [ $nproc == 512 ]; then
                    cp ~/run.sh . 
                    cp ~/rr     . 
                    sed -i '/case/c\case=b'$ele'' run.sh
                    sed -i '/qsub/c\ qsub -n 8 ./run.sh 64' rr
                echo 'hello',$proc,$nproc,$ele
               fi
	    fi
	fi
	done
    done
done
