/**
 * This file contains common functions that are needed in MPI-IO schemas, mainly
 * get file name, start/end timing and write io trace.
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "mpiio_util.h"

char nmFilename[kMaxPathLen];
//char filename[kMaxPathLen];
char mFilename[kMaxPathLen];
char rstFilename[kMaxPathLen];
char rbFilename[kMaxPathLen];
char rbasciiFilename[kMaxPathLen];
char rbnmmFilename[kMaxPathLen];

char thefilename[kMaxPathLen]; // keep filename of current io_option

char path[kMaxPathLen];

double start_time, end_time;
double overall_time = 0;
double io_time = 0;
double file_io_time = 0;

int trace_ioop = -1;
int trace_nf = -1;

void getfilename_(int *id, int *nid, int io_option)
{
  double start_getfilename = MPI_Wtime();
  if(DEBUG_FLAG) printf("io_option = %d, numGroups = %d\n", io_option, numGroups);

	DIR* dir = NULL;
	char ext0[100];
	char ext1[100];

	/*printf( "\n  nid:: %d\n", *nid);*/
	memset((void*)filename, 0, kMaxPathLen);
	memset((void*)rstFilename, 0, kMaxPathLen);
	memset((void*)mFilename, 0, kMaxPathLen);
	memset((void*)rbFilename, 0, kMaxPathLen);
	memset((void*)rbnmmFilename, 0, kMaxPathLen);
	memset((void*)nmFilename, 0, kMaxPathLen);
	memset((void*)thefilename, 0, kMaxPathLen);
	//char path[kMaxPathLen];
	memset((void*)path, 0, kMaxPathLen);

	sprintf(path, kOutputPath);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &mysize);

	// if it's local, the "vtk" dir is already created,
	// simply put everything in this dir (mostly for debug)
	if(strcmp(kOutputPath, kStrLocal) == 0) {
	//if(0) {
		if(myrank == 1) printf("Output files will be in local dir %s\n", path);
		sprintf(filename, "%s/binary-NN-p%.6d-t%.5d.vtk", path, *nid, *id);
		sprintf(rstFilename, "%s/restart-mpi-binary-N1-t%.5d.vtk",path, *id);
		sprintf(mFilename, "%s/mpi-binary-N1-t%.5d.vtk",path, *id);
		sprintf(rbFilename, "%s/mpi-binary-NM1-t%.5d.vtk", path, *id);
		sprintf(rbasciiFilename, "%s/mpi-ascii-NM1-t%.5d.vtk", path, *id);
		sprintf(nmFilename, "%s/mpi-binary-NM-p%.6d-t%.5d.vtk",
						path, groupRank, *id);
		if(io_option == 8) 
      sprintf(rbnmmFilename, "%s/mpi-binary-NMM-p%.6d-t%.5d.vtk", 
              path, groupRank, *id);
    else if (io_option == 18) 
      sprintf(rbnmmFilename, "%s/mpi-binary-NMM-thread-p%.6d-t%.5d.vtk", 
              path, groupRank, *id);

	}
	else if (kOutputPath != NULL) {
		//rank 0 create top level dir
		if(myrank == 0) {
      // create top-level dir if not exist
			dir = opendir(path);
			//if non-exist, create it
			if(dir == NULL) {
        printf("\nOutput path %s does not exist, create new one\n");
				int status = mkdir(path, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
				if(status != 0) {
					printf("can't create dir for top-level output dir at %s\n", path);
					exit(1);
				}
			}
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
			//for io_option 5 ,8 and 18, it have a NM layer dir
			if(io_option == 5 || io_option == 8 || io_option == 18) {
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
			if(DEBUG_FLAG) printf("output path is %s\n", path);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(path, sizeof(path), MPI_CHAR, 0, MPI_COMM_WORLD);

		if(io_option ==99) {
			//printf("io-option");
			sprintf(rstFilename, "%s/restart-%d-proc-mpi-binary-N1-t%.5d.vtk",
						path, mysize, *id);
		}
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
		else if(io_option == -6) {
		sprintf(rbasciiFilename, "%s/%d-proc-mpi-ascii-rbIO-NM1-t%.5d.vtk",
						path, mysize,*id);
		}
		//generating NM files, create dir for them
		else if(io_option == 5 || io_option == 8 || io_option == 18) {
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
      else if(io_option == 18) {
				sprintf(rbnmmFilename,
						"%s/%d-proc-mpi-binary-rbIO-NMM-thread-p%.6d-t%.5d.vtk",
						path,  mysize, groupRank, *id);
			}
			else if(io_option == 5) {
				sprintf(nmFilename,
								"%s/%d-proc-mpi-binary-coIO-NM-p%.6d-t%.5d.vtk",
						path,  mysize,groupRank, *id);
			}
		}// end of if 5 or 8 or 18
	}
	else {
		printf("error: the kOutputPath %s does not match anything in getfilename(), please check again\n", kOutputPath);
		exit(5);
	}
	adjust_endian();

  double end_getfilename = MPI_Wtime();
  if(myrank == 0) 
    printf("rank 0 getfilename (and create dir if non-exist) takes %lf sec.\n", 
           end_getfilename - start_getfilename);
}

#ifdef UPCASE
void STARTTIMING()
#elif  IBM
void starttiming()
#else
void starttiming_()
#endif
{
  start_time =  MPI_Wtime();
}

#ifdef UPCASE
void ENDTIMING()
#elif  IBM
void endtiming()
#else
void endtiming_()
#endif
{
  end_time = MPI_Wtime();
  overall_time = end_time - start_time;
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

	char tracefname[128];
	int formatparam = *fparam;
	int iostep = *piostep;

	memset((void*)tracefname, 0, 128);
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
  double io_time_max = 0.0;
  double file_io_max = 0.0;
	if( formatparam == 2 || formatparam == 3 || formatparam == 4 || formatparam == 5)
	{
		MPI_Comm_size(MPI_COMM_WORLD, &mysize);
		MPI_Allreduce(  &overall_time, &overall_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
		MPI_Allreduce(  &overall_time, &overall_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
		MPI_Allreduce(  &overall_time, &overall_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		overall_avg = overall_sum / mysize;
	}
	else if(formatparam == 6 || formatparam == -6 || formatparam == 8 || formatparam == 18)
	{
		if(mySpecies == 1)
		{
		MPI_Allreduce(  &overall_time, &overall_min, 1, MPI_DOUBLE, MPI_MIN, localcomm);
		MPI_Allreduce(  &overall_time, &overall_max, 1, MPI_DOUBLE, MPI_MAX, localcomm);
		MPI_Allreduce(  &overall_time, &overall_sum, 1, MPI_DOUBLE, MPI_SUM, localcomm);
		overall_avg = overall_sum / localsize;

		MPI_Allreduce(  &file_io_time, &file_io_max, 1, MPI_DOUBLE, MPI_MAX, localcomm);
		MPI_Allreduce(  &io_time, &io_time_max, 1, MPI_DOUBLE, MPI_MAX, localcomm);
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
		
    printf("**************************************\n");
		printf("I/O time (io_step=%d) stats: overall avg = %lf sec, min = %lf sec, max = %lf sec "
           "(io_max = %lf sec, file_io_max = %lf sec, wtick=%lf sec),"
					 "checkpoint file path is %s, machine is %s, io_option = %d, num_groups = %d "
           "(DEBUG_FLAG=%d, COMPUTE_TRACE_FLAG=%d).\n",
					 io_step, overall_avg, overall_min, overall_max, io_time_max, file_io_max, MPI_Wtick(),
           path, mach_name, formatparam, numGroups,
           DEBUG_FLAG, COMPUTE_TRACE_FLAG);
    printf("**************************************\n");
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if(IOTRACE_FLAG != 1)
		return;

  // write the actual file
  if (0) {
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
#ifdef UPCASE
void PASS_IO_PARAMS(int *param1, int* param2)
#elif  IBM
void pass_io_params(int *param1, int* param2)
#else
void pass_io_params_(int *param1, int* param2)
#endif
{
  trace_ioop= *param1;
  trace_nf = *param2;
  printf("io_option = %d, nfiles = %d", trace_ioop, trace_nf);
}

#ifdef UPCASE
void WRITECOMPUTETRACE(int* pcompstep, double* pdtime, double* pcpu_t)
#elif  IBM
void writecomputetrace(int* pcompstep, double* pdtime, double* pcpu_t)
#else
void writecomputetrace_(int* pcompstep, double* pdtime, double* pcpu_t)
#endif
{
	//printf("format param is %d, iostep is %d, dtime = %lf\n", (int)*fparam, *pcompstep, *pcpu_t);
	if(COMPUTE_TRACE_FLAG != 1)
		return;

	char tracefname[kMaxPathLen];
	int formatparam = trace_ioop;
  int nfile = trace_nf;
	int stepnum = *pcompstep;
  double dtime = *pdtime;
  double cpu_t = *pcpu_t;

	int temp_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &temp_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &mysize);

	MPI_Barrier(MPI_COMM_WORLD);

	memset((void*)tracefname, 0, kMaxPathLen);
	//sprintf(tracefname, "%s/compute-trace-%d-proc-ioop-%d-nf-%d-t%.5d.dat", 
   //       kOutputPath, mysize, formatparam, nfile, stepnum);

  sprintf(tracefname, "%s/../compute-trace-t%.5d.dat", path, stepnum);

  //if(myrank == 0) printf("compute filename = %s", tracefname);
  // write the actual file
  if (1) {
		MPI_File timefile;
		int rc;
		rc = MPI_File_open(MPI_COMM_WORLD, tracefname,
									 		MPI_MODE_CREATE | MPI_MODE_WRONLY , MPI_INFO_NULL, &timefile);
    if(rc) {
      if(temp_rank == 0) printf("Unble to open file %s! \n", tracefname);
    }

		char mytime[128];
    memset(mytime, 0, 128);
		sprintf(mytime, "%10d %10.3lf %10.3lf\n",
						temp_rank, dtime, cpu_t);

    int len = strlen(mytime); 
    //printf("str len = %d\n", len);

		long long offsets = temp_rank * len ;
		MPI_Status write_data_status;

		MPI_File_write_at_all_begin(timefile,
													 			offsets,
																mytime,
																len,
																MPI_CHAR);
		MPI_File_write_at_all_end(timefile,
															mytime,
															&write_data_status);
		MPI_File_close( & timefile );
	}
}
