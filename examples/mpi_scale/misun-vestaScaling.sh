#!/bin/bash

# this is for MPI runs on Vesta/Cetus/Mira

N=7

for proc in 'MPI' ; do
    for nproc in 1 2 4 8 16 32 64 128; do
	for ele in 1 2 4 8 16 32 64 128 256 512 1024 2048; do
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
	    	cp ../mpi_scale/data/b$ele.* .
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
            ../../bin/makenekmpi

	    #Submit job 
	    if    [ $ele -lt 32 ]; then
	       if    [ $nproc -le 32 ]; then
 		../../bin/nek box $nproc 1 
		echo 'hello',$proc,$nproc,$ele
	       elif  [ $nproc == 64 ]; then
 		../../bin/nek box $nproc 2 
		echo 'hello',$proc,$nproc,$ele
	       elif  [ $nproc == 128 ]; then
 		../../bin/nek box $nproc 4 
		echo 'hello',$proc,$nproc,$ele
               fi
            else   
	       if    [ $nproc -le 32 ]; then
 		../../bin/nek b$ele $nproc 1 
		echo 'hello',$proc,$nproc,$ele
	       elif  [ $nproc == 64 ]; then
 		../../bin/nek b$ele $nproc 2 
		echo 'hello',$proc,$nproc,$ele
	       elif [ $nproc == 128 ]; then
 		../../bin/nek b$ele $nproc 4 
		echo 'hello',$proc,$nproc,$ele
	       fi
	    fi
	fi
	done
    done
done
