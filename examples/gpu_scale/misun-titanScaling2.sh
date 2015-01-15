#!/bin/bash

# this is for MPI runs on Vesta/Cetus/Mira

N=7

for proc in 'MPI'; do
    for nproc in 1 2 4 8 16 32 64 128; do
	for ele in 1 2 4 8 16 32 64 128 256 512 1024 2048 4096; do
            if [ $ele -ge $nproc ]; then
            
            
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
	    	cp ../gpu_scale/nekcem .
    	        cp ../gpu_scale/box.* .
                
	    elif [ $proc == 'MPI' ]; then
	    	cp ../mpi_scale/nekcem .
	    	cp ../mpi_scale/box.* .
	    else
	    	echo 'Must be MPI or GPU '
	    	exit
	    fi
	    
	    #Change box file
	    sed -i '/(number of elements/c\-1 -1 -'$ele' (number of elements in x,y,z)' box.box
            sed -i '/NEL,NDIM,NELV/c\   '$ele'    -3   '$ele'           NEL,NDIM,NELV' box.rea

	    #Submit job 
            if [ $proc == 'GPU' ]; then
                if [ $nproc -gt 8 ]; then
                    break
                fi
#                module swap PrgEnv-pgi PrgEnv-cray
#                module load craype-accel-nvidia35
#                module load perftools
#                ../../bin/cleanall   
#                ../../bin/makenekgpu -a titan-cce-acc
                ../../bin/nekgpu box $nproc $nproc
                echo "Hello",$nproc,$nproc,$ele
            elif [ $proc == 'MPI' ]; then 
#                module swap PrgEnv-pgi PrgEnv-cray
#                module unload craype-accel-nvidia35
#                module unload perftools
#                ../../bin/cleanall   
#                ../../bin/makenekmpi -a titan-cce-mpi
#Calculate number of nodes needed
                nnode=$((($nproc-1)/16 + 1))
                #Number of processors followed by number of nodes
                ../../bin/nek box $nproc $nnode 
                echo "Hello",$nnode,$nproc,$ele
            else
                echo 'Must be MPI or GPU'
               exit
            fi
            fi
	done
    done
done
