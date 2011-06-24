#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "mpiio_util.h"

char filename[100];
char mFilename[100];
char rbFilename[100];
char rbnmmFilename[128];
char nmFilename[128];

long long start_time, end_time;
double overall_time;

void getfilename_(int *id, int *nid, int io_option)
{
	DIR* dir = NULL;
	char ext0[100];
	char ext1[100];

	/*printf( "\n  nid:: %d\n", *nid);*/
	memset((void*)filename, 0, 100);
	memset((void*)mFilename, 0, 100);
	memset((void*)rbFilename, 0, 100);
	memset((void*)rbnmmFilename, 0, 128);
	memset((void*)nmFilename, 0, 128);
	char path[128];
	memset((void*)path, 0, 128);

	sprintf(path, kOutputPath);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &mysize);

	// if it's local, the "vtk" dir is already created,
	// simply put everything in this dir (mostly for debug)
	if(strcmp(kOutputPath, kStrLocal) == 0) {
		if(myrank == 1) printf("Output files will be in local dir %s\n", path);
		sprintf(filename, "%s/binary-NN-p%.6d-t%.5d.vtk", path, *nid, *id);
		sprintf(mFilename, "%s/mpi-binary-N1-t%.5d.vtk",path, *id);
		sprintf(rbFilename, "%s/mpi-binary-NM1-t%.5d.vtk", path, *id);
		sprintf(rbasciiFilename, "%s/mpi-ascii-NM1-t%.5d.vtk", path, *id);
		sprintf(rbnmmFilename, "%s/mpi-binary-NMM-p%.6d-t%.5d.vtk",
						path, groupRank, *id);
		sprintf(nmFilename, "%s/mpi-binary-NM-p%.6d-t%.5d.vtk",
						path, groupRank, *id);
	}
	else if(strcmp(kOutputPath, kStrFs0Misun) == 0 ||
			strcmp(kOutputPath, kStrFs0Fuji) == 0) {
		//rank 0 create top level dir
		if(myrank == 0) {
			//create NP/IO_OPTION directory first
			sprintf(path, "%s/%d", path, mysize);
			dir = opendir(path);
			//if non-exist, create it
			if(dir == NULL) {
				int status = mkdir(path, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
				if(status != 0) {
					printf("can't create dir for /NP\n");
					exit(1);
				}
			}
			else {
				assert(closedir(dir) == 0);
			}
			//now it exists, let's add /IO_OPTION to it
			sprintf(path, "%s/%d", path, io_option);
			dir = opendir(path);
			//if non-exist, create it
			if(dir == NULL) {
				int status = mkdir(path, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
				if(status != 0) {
					printf("can't create dir for /io_option\n");
					exit(2);
				}
			}
			else {
				assert(closedir(dir) == 0);
			}
			//for io_option 5 and 8, it have a NM layer dir
			if(io_option == 5 || io_option == 8) {
				sprintf(path, "%s/%d", path, numGroups);
				dir = opendir(path);
				if(dir == NULL) {
					int status = mkdir(path, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
					if(status != 0) {
						printf("can't create dir for /NM\n");
						exit(3);
					}
				}
				else {
					assert(closedir(dir) == 0);
				}
			}
			printf("output path is %s\n", path);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(path, sizeof(path), MPI_CHAR, 0, MPI_COMM_WORLD);

		if(io_option == 3) {
			sprintf(filename, "%s/%d-proc-binary-NN-p%.6d-t%.5d.vtk",
						path, mysize, *nid, *id);
		}
		if(io_option == 4) {
			sprintf(mFilename, "%s/%d-proc-mpi-binary-N1-t%.5d.vtk",
						path, mysize, *id);
		}
		//for this case, only 1 file generated, so no seperate dir
		else if(io_option == 6) {
			sprintf(rbFilename, "%s/%d-proc-mpi-binary-rbIO-NM1-t%.5d.vtk",
							path, mysize, *id);
		}
		//ascii version of io_option = 6, same config
		else if(io_option == 7) {
		sprintf(rbasciiFilename, "%s/%d-proc-mpi-ascii-rbIO-NM1-t%.5d.vtk",
						path, mysize,*id);
		}
		//generating NM files, create dir for them
		else if(io_option == 5 || io_option == 8) {
			sprintf(path, "%s/%d", path, groupRank);
			if(mySpecies == 1) {
				dir = opendir(path);
				if(dir == NULL) {
					int status = mkdir(path, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
					if(status != 0) {
						printf("can't create dir for each group\n");
						exit(4);
					}
				}
				else {
					assert(closedir(dir) == 0);
				}
			}
			if(io_option == 8) {
				sprintf(rbnmmFilename,
						"%s/%d-proc-mpi-binary-rbIO-NMM-p%.6d-t%.5d.vtk",
						path,  mysize, groupRank, *id);
			}
			else if(io_option == 5) {
				sprintf(nmFilename,
								"%s/%d-proc-mpi-binary-syncIO-NM-p%.6d-t%.5d.vtk",
						path,  mysize,groupRank, *id);
			}
		}// end of if 5 or 8
	}
	else {
		printf("error: the kOutputPath does not match anything in getfilename(), please check again\n");
		exit(5);
	}
	adjust_endian();
}

#ifdef UPCASE
void STARTTIMING()
#elif  IBM
void starttiming()
#else
void starttiming_()
#endif
{
        start_time = rdtsc();
}

#ifdef UPCASE
void ENDTIMING()
#elif  IBM
void endtiming()
#else
void endtiming_()
#endif
{
        end_time = rdtsc();
	overall_time = (double) (end_time - start_time)/ (BGP_FREQ) ;
        if(IOTIMER_FLAG)
	{
//		if(myrank == 0)
//			printf("\noverall I/O time is %lf seconds \n", overall_time);
	}
}

#ifdef UPCASE
void WRITEIOTRACE(int *fparam, int* piostep)
#elif  IBM
void writeiotrace(int *fparam, int* piostep)
#else
void writeiotrace_(int *fparam, int* piostep)
#endif
{
	//printf("format param is %d, iostep is %d\n", (int)*fparam, *piostep);
	if(IOTRACE_FLAG != 1)
		return;

	char tracefname[128];
	int formatparam = *fparam;
	int iostep = *piostep;

	memset((void*)tracefname, 128, '\0');
/*
	if(formatparam == 2) sprintf(tracefname, "ascii-NN-iotrace");
	else if(formatparam == 3) sprintf(tracefname, "binary-NN-iotrace");
	else if(formatparam == 4) sprintf(tracefname, "mpi-binary-N1-iotrace");
	else if(formatparam == 6) sprintf(tracefname, "mpi-binary-NM1-iotrace");
	else if(formatparam == -6) sprintf(tracefname, "mpi-ascii-NM1-iotrace");
	else if(formatparam == 8) sprintf(tracefname, "mpi-binary-NMM-iotrace");
	sprintf(tracefname, "%s-t%.5d.dat", tracefname, iostep);
*/
	sprintf(tracefname, "iotrace-t%.5d.dat", iostep);

	double overall_max, overall_min, overall_avg, overall_sum;
	if( formatparam == 2 || formatparam == 3 || formatparam == 4 || formatparam == 5)
	{
		MPI_Comm_size(MPI_COMM_WORLD, &mysize);
		MPI_Allreduce(  &overall_time, &overall_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
		MPI_Allreduce(  &overall_time, &overall_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
		MPI_Allreduce(  &overall_time, &overall_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		overall_avg = overall_sum / mysize;
	}
	else if(formatparam == 6 || formatparam == -6 || formatparam == 8)
	{
		if(mySpecies == 1)
		{
		MPI_Allreduce(  &overall_time, &overall_min, 1, MPI_DOUBLE, MPI_MIN, localcomm);
		MPI_Allreduce(  &overall_time, &overall_max, 1, MPI_DOUBLE, MPI_MAX, localcomm);
		MPI_Allreduce(  &overall_time, &overall_sum, 1, MPI_DOUBLE, MPI_SUM, localcomm);
		overall_avg = overall_sum / localsize;

		}
		else if(mySpecies == 2)
		{
			overall_time = 0;
			overall_min = 0;
			overall_max = 0;
			overall_avg = 0;
		}
	}

	int temp_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &temp_rank);

	if(temp_rank == 0) {
		printf("I/O time - avg = %lf seconds, max = %lf seconds ,"
					 "restart file dir is %s(show fs0 or local)\n",
					 overall_avg, overall_max, filename);
	}
	MPI_Barrier(MPI_COMM_WORLD);
        if (0)
	{
		MPI_File timefile;
		int rc;
		rc = MPI_File_open(MPI_COMM_WORLD, tracefname,
									 		MPI_MODE_CREATE | MPI_MODE_WRONLY , MPI_INFO_NULL, &timefile);

		char mytime[128];
		sprintf(mytime, "\n%10d %10.3lf %10.3lf %10.3lf %10.3lf ",
						temp_rank, overall_time, overall_avg, overall_min, overall_max);

		long long offsets = temp_rank * 56 ;
		MPI_Status write_data_status;

		MPI_File_write_at_all_begin(timefile,
													 			offsets,
																mytime,
																56,
																MPI_CHAR);
		MPI_File_write_at_all_end(timefile,
															mytime,
															&write_data_status);
		MPI_File_close( & timefile );
	}
}



