!/bin/bash

# this is for MPI runs on Vesta/Cetus/Mira/Theta

for proc in 'MPI' ; do
    for nproc in 1 2 4 8 16 32 64 128 256 512; do
        for ele in 1 2 4 8 16 32 64 128 256 512 1024 2048 4096; do
            dir=${proc}_${nproc}_${ele}

            #Make new directory
            if [ "$ele" -ge "$nproc" ]; then

            cp ${dir}/*output ./mpi_scale/data/smallOutput/theta/MPI/N14/theta_${proc}_${nproc}_${ele}
            fi
        done
    done
done

