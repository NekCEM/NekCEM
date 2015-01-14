#!/bin/bash

# this is for MPI runs on Vesta/Cetus/Mira

N=7

for proc in 'MPI' 'GPU'; do
    for nproc in 1 2 4 8 16 32 64 128; do
	for ele in 1 2 4 8 16 32 64 128 256 512 1024; do
	    dir=${proc}_${nproc}_${ele}

	    #Make new directory
	    mkdir ../${dir}
	    cd ../${dir}
	    echo $proc

	    #copy in the right files
            if [ $proc == 'GPU' ]; then
                if [ $nproc -gt 8 ]; then
                    break
                fi
	    	cp ../gpu_scale/SIZEu .
	    	cp ../gpu_scale/box.* .

	    elif [ $proc == 'MPI' ]; then
	    	cp ../mpi_scale/SIZEu .
	    	cp ../mpi_scale/box.* .
	    else
	    	echo 'Must be MPI or GPU '
	    	exit
	    fi
	    
	    #Change box file
	    sed -i '/(number of elements/c\-1 -1 -'$ele' (number of elements in x,y,z)' box.box
            sed -i '/polynomial order/c\      parameter (lxi =   '$N')  ! polynomial order' SIZEu
            sed -i '/number of mpi/c\      parameter (lp  =   '$nproc')  ! number of mpi' SIZEu
            sed -i '/number of lelg/c\      parameter (lelg=   '$ele')  ! number of lelg' SIZEu
            sed -i '/number of lelz/c\      parameter (lelx =1,lely =lelx ,lelz =  '$ele')  ! number of lelz' SIZEu
            sed -i '/NEL,NDIM,NELV/c\   '$ele'    -3   '$ele'           NEL,NDIM,NELV' box.rea

	    #Submit job 
            if [ $proc == 'GPU' ]; then
                if [ $nproc -gt 8 ]; then
                    break
                fi
                ../../bin/cleanall   
                ../../bin/makenekgpu
                ../../bin/nekgpu box $nproc $nproc
 
            elif [ $proc == 'MPI' ]; then 
                ../../bin/cleanall   
                ../../bin/makenekmpi
                ../../bin/nek box 1 $nproc
            else
                echo 'Must be MPI or GPU'
               exit
            fi
	done
    done
done
