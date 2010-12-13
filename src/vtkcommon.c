#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "vtkcommon.h"

/*#ifdef NEED_TRAILING_UNDERSCORE
   #define FORTRAN(SUBROUTINE_NAME) SUBROUTINE_NAME##_
#else
  #define FORTRAN(SUBROUTINE_NAME) SUBROUTINE_NAME
#endif
FORTRAN(getfieldname) */

/*
#ifdef NEED_TRAILING_UNDERSCORE
  #define FORTRAN(SUBROUTINE_NAME) SUBROUTINE_NAME##_
#else
  #define FORTRAN(SUBROUTINE_NAME) SUBROUTINE_NAME
#endif
*/

/*
FORTRAN(openfile)
FORTRAN(closefile)
FORTRAN(writeheader)  
FORTRAN(writenodes)   
FORTRAN(write2dcells) 
FORTRAN(write3dcells) 
FORTRAN(writefield)   
*/

int Little_endian = -1; 

char filename[100];
char mFilename[100];
char rbFilename[100];
char rbnmmFilename[128];
char nmFilename[128];

long long start_time, end_time;
double overall_time;

int DEBUG_FLAG = 0;
int IOTIMER_FLAG = 1;
int IOTRACE_FLAG = 1;


/**************************/


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

void getfieldname_( int i, char *name )
{
   int id = i;
   switch( id )
   {
	   case 1:
		   strcpy(name, "H ");
		   break;
	   case 2:
		   strcpy(name, "E ");
		   break;
	   case 11:
		   strcpy(name, "errH");
		   break;
	   case 12:
		   strcpy(name, "errE");
		   break;
	   case 13:
		   strcpy(name, "solH");
		   break;
	   case 14:
		   strcpy(name, "solE");
		   break;
	   case 21:
		   strcpy(name, "totH");
		   break;
	   case 22:
		   strcpy(name, "totE");
		   break;
	   case 23:
		   strcpy(name, "scatH");
		   break;
	   case 24:
		   strcpy(name, "scatE");
		   break;
	   case 25:
		   strcpy(name, "incH");
		   break;
	   case 26:
		   strcpy(name, "incE");
		   break;
	   case 27:
		   strcpy(name, "timeavgE");
		   break;
	   case 28:
		   strcpy(name, "timeavgIE");
		   break;
	   case 29:
		   strcpy(name, "poyVec");
		   break;
	   case 30:
		   strcpy(name, "poytVec");
		   break;
	   case 31:
		   strcpy(name, "poyVecI");
		   break;
	   case 32:
		   strcpy(name, "avgVec");
		   break;
	   case 33:
		   strcpy(name, "avgtVec");
		   break;
	   case 34:
		   strcpy(name, "avgVecI");
		   break;
	   case 201:
		   strcpy(name, "region");
		   break;
	   case 202:
		   strcpy(name, "face");
		   break;
	   case 203:
		   strcpy(name, "No ");
		   break;
	   case 204:
		   strcpy(name, "PEC");
		   break;
	   case 205:
		   strcpy(name, "PML");
		   break;
	   case 206:
		   strcpy(name, "P  ");
		   break;
	   case 101:
		   strcpy(name, "potE");
		   break;
	   case 102:
		   strcpy(name, "spotE");
		   break;
   }
}

//this function detects machine's endianness
void adjust_endian()
{
	int endian_int = 1;
	char* pchar = &endian_int;
	//	for(int i = 0 ; i < 4; i++) printf(" -%d ",(int) *(pchar+i));
	if(* (pchar+3)  == 1) Little_endian = 0;
	else Little_endian = 1;

	if(DEBUG_FLAG == 3)
	{
		if(Little_endian == 0)printf("I'm big endian\n");
		else if(Little_endian == 1) printf("I'm little endian\n");
		else printf("Endianness Error!!!\n"); 
	}
}

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
	MPI_Comm_size(MPI_COMM_WORLD, &mysize);

	if(strcmp(kOutputPath, kStrLocal) == 0) {
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
					printf("can't create dir for /NP");
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
					printf("can't create dir for /io_option");
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
						printf("can't create dir for /NM");
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
					close(dir);
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
	
	adjust_endian();
}


int swap_int_byte(int *n)
{
	//if it's little endian, don't have to do extra byte swapping
	if(Little_endian == 1)
	{
		unsigned char *cptr,tmp;

		cptr = (unsigned char *)n;
		tmp = cptr[0];
		cptr[0] = cptr[3];
		cptr[3] = tmp;
		tmp = cptr[1];
		cptr[1] = cptr[2];
		cptr[2] = tmp;
		//printf("I'm little endian\n");
	}
  return 0;
}

int swap_float_byte(float *n)
{
if(Little_endian == 1)
{
  unsigned char *cptr,tmp;

  cptr = (unsigned char *)n;
  tmp  = cptr[0];
  cptr[0] = cptr[3];
  cptr[3] = tmp    ;
  tmp     = cptr[1];
  cptr[1] = cptr[2];
  cptr[2] = tmp    ;
}
  return 0;
}

//added by Jing Fu at 2010-7-22
int swap_long_long_byte(long long *n)
{
	if(Little_endian == 1)
	{
		unsigned char *cptr, tmp;
		cptr = (unsigned char*)n;
		tmp = cptr[0];
		cptr[0] = cptr[7];
		cptr[7] = tmp;
		
		tmp = cptr[1];
		cptr[1] = cptr[6];
		cptr[6] = tmp;

		tmp = cptr[2];
		cptr[2] = cptr[5];
		cptr[5] = tmp;
		
		tmp = cptr[3];
		cptr[3] = cptr[4];
		cptr[4] = tmp;
	}
	return 0;
}



