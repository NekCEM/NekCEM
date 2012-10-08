/**
 * This file contains functions for rbIO (normal and thread version).
 *
 * It requires MPI.
 */
#include <stdio.h>
#include <pthread.h>
#include <errno.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <stdint.h>

#include "mpiio_util.h"

#define HYBRID_PTHREAD
int THREAD = 0;

extern MPI_File mfile;
int64_t testint64;
int fieldSizeLimit ;

//char rbFilename[100];
//char rbasciiFilename[100];
//char mfBuffer[ 4* ONE_MILLION];
//long long mfBufferCur  = 0;
//long long mfileCur     = 0;
//long long fieldSizeSum = 0;

char *sendBuffer;
int   sendBufferCur = 0;

/*global rank and size*/
int myrank, mysize           ;
int localrank, localsize     ; //rank and size in localcomm
MPI_Comm localcomm, groupcomm; //worker comm and write comm
int numGroups, groupSize, rankInGroup, groupRank,  mySpecies, numFields, ifield;

file_t*		file = NULL;

char*			writerBuffer ;
char**		recvBuffers  ;
char*			recvmsgBuffer;
int*			iSize;
long long writerBufferCur = 0 ;
int       writerBufferSize    ;
int       recvmsgBufferCur = 0;

int first_init   = 0;
int AUGMENT_FLAG = 0;
int ASCII_FLAG   = 0;
int io_option		 = 0;

int INT_DIGITS   = 10;
int FLOAT_DIGITS = 18;
int LONG_LONG_DIGITS = 18;

double io_time_start, io_time_end;

int i, j; //so that intel c compiler won't complain def in a loop

void set_io_option( int option)
{
  if (option == 5) {
    ASCII_FLAG = 3;
    io_option = 5;
  }
  else if (option == 6) {
    ASCII_FLAG = 0;
    io_option = 6;
  }
  else if (option == -6) {
    ASCII_FLAG = 1;
    io_option = 7;
  }
  else if (option == 8) {
    ASCII_FLAG = 2;
    io_option = 8;
  }
  else if (option == 18) {
    ASCII_FLAG = 2;
    io_option = 18;
    THREAD = 1; // threaded rbIO NMM
  }
  else {
    printf("ERROR: wrong io_option value passed to set_io_option!\n");
    exit(1);
  }
}

void set_ascii_true ()
{
	ASCII_FLAG = 1;
	if (DEBUG_FLAG) printf("setting ascii flag to true\n");
}

void set_ascii_nmm ()
{
	ASCII_FLAG = 2;
	if (DEBUG_FLAG) printf("setting ascii flag to NMM\n");
}

void set_ascii_nm ()
{
	ASCII_FLAG = 3;
	if (DEBUG_FLAG) printf("setting ascii flag to NM\n");
}

#ifdef UPCASE
void FREE_RBIO_BUFFER ()
#elif  IBM
void free_rbio_buffer ()
#else
void free_rbio_buffer_()
#endif
{
	free(mfBuffer);
	free(sendBuffer);
	for( i = 0; i < groupSize; i++)
		free (recvBuffers[i]) ;
	free(recvBuffers);
	free(iSize);
	free(recvmsgBuffer);
	free(writerBuffer);
}

void smartCheckGroupSize(int *numgroups)
{
//	int IDEAL_SIZE = 64;
//	int SIZE_UPPER_BOUND = 64;
	if( (mysize/(*numgroups)) > GROUP_SIZE_UPPER_BOUND)
	{
		*numgroups = mysize/GROUP_SIZE_IDEAL;
		if(myrank == 0)printf("changed numGroups to %d\n", *numgroups);
	}
}

#ifdef UPCASE
void INITRBIO (int *numgroups, int* numfields, int* maxnumnodes)
#elif  IBM
void initrbio (int *numgroups, int* numfields, int* maxnumnodes)
#else
void initrbio_(int *numgroups, int* numfields, int* maxnumnodes)
#endif
{
	// this if is only executed once - first time
	if (first_init == 0) {
		first_init = 1; //only init for first time

		MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
		MPI_Comm_size(MPI_COMM_WORLD, &mysize);

		smartCheckGroupSize(numgroups);

		numGroups = *numgroups    ; //printf("numgroups is %d*****\n", *numgroups);
		numFields = *numfields + 3; //+3 due to node coord, cell numbering and celltype

		groupSize   = mysize / numGroups;
		groupRank   = myrank / groupSize;
		rankInGroup = myrank % groupSize;

		/*
		 *upper bound of single field size(include nodes info and 3d cells)
		 *for ascii case, float has 18 digits * 3, int has 7 digits * 10, so take int size
		 */
		if      (ASCII_FLAG == 0 || ASCII_FLAG == 2 || ASCII_FLAG == 3)
			fieldSizeLimit = 10 *  sizeof(int) * (*maxnumnodes) + 1024;
		else if (ASCII_FLAG == 1)
			fieldSizeLimit = 10 * (INT_DIGITS + 1) * (*maxnumnodes) + 1024;

		//ASCII_FLAG=1 is the real ASCII case (NM1, i.e. rbIO 1 file)
		if      (ASCII_FLAG == 0 || ASCII_FLAG == 1 || ASCII_FLAG == 3)
			writerBufferSize = groupSize * fieldSizeLimit;
		else if (ASCII_FLAG == 2)
			writerBufferSize = WRITERBUFFERSIZE;

		mfBuffer   = (char*) malloc(sizeof(char) * fieldSizeLimit);
		sendBuffer = (char*) malloc(sizeof(char) * fieldSizeLimit);

		// if this proc is a writer in its sub-group
		// writer species is 1, worker is 2
		if (rankInGroup == 0)
		{
			mySpecies  = 1;
			recvBuffers= (char**)malloc(sizeof(char*) * groupSize);
			for( i = 0; i < groupSize; i++)
				recvBuffers[i] = (char*) malloc(sizeof(char) * fieldSizeLimit);
			iSize = (int*) malloc(sizeof(int) * groupSize);
			recvmsgBuffer = (char*) malloc(sizeof(char) * fieldSizeLimit  );
			writerBuffer  = (char*) malloc(sizeof(char) * writerBufferSize);
			if (writerBuffer == NULL)
			{
				printf("not enough memory for writerBuffer - quitting!!\n");
				exit(1);
			}
			init_file_struc();
		}
		else mySpecies = 2;

		MPI_Comm_split(MPI_COMM_WORLD, mySpecies, myrank, &localcomm);
		MPI_Comm_split(MPI_COMM_WORLD, groupRank, myrank, &groupcomm);
		MPI_Comm_rank(localcomm, &localrank);
		MPI_Comm_size(localcomm, &localsize);

		if (DEBUG_FLAG) printf("myrank is %d, rankInGroup = %d, localrank is %d, numGroups is %d, fieldSizeLimit is %d, maxnumnodes is %d\n", myrank, rankInGroup, localrank, numGroups, fieldSizeLimit, *maxnumnodes);
	}// end of if first time

	// if it's writer, reset file struc for every write
	if(mySpecies == 1) reset_file_struc();

	ifield      = 0;
	mfBufferCur = 0;
	mfileCur    = 0; //reset file position for new file
}

#ifdef UPCASE
void OPENFILE6(  int *id, int *nid)
#elif  IBM
void openfile6(  int *id, int *nid)
#else
void openfile6_(  int *id, int *nid)
#endif
{
  io_time_start = MPI_Wtime();

	getfilename_(id, nid, io_option);

	if(ASCII_FLAG == 3) // i.e. io_option == 5, sync M files
	{
		int rc = MPI_File_open(groupcomm, nmFilename, MPI_MODE_CREATE | MPI_MODE_RDWR , MPI_INFO_NULL, &mfile);
		if(rc){
			printf("Unable to create unique file %s in openfile in rank %d\n", nmFilename, myrank);
			fflush(stdout);
		}

	}
	else if(mySpecies == 1)
	{
		/* parallel here*/

		int rc;
		if(ASCII_FLAG == 0) // i.e. io_option = 6, rbIO 1 file
		{
			rc = MPI_File_open(localcomm, rbFilename, MPI_MODE_CREATE | MPI_MODE_RDWR , MPI_INFO_NULL, &mfile);
           if (myrank ==0)  printf("openfile6: rc ---2 %d\n,rc");
			if(rc){
				printf("Unable to create shared file %s in openfile\n", rbFilename);
				fflush(stdout);
			}
		}
		else if(ASCII_FLAG == 1) // i.e. io_option = 7, true ASCII rb 1 file
		{
			rc = MPI_File_open(localcomm, rbasciiFilename, MPI_MODE_CREATE | MPI_MODE_RDWR , MPI_INFO_NULL, &mfile);
			if(rc){
				printf("Unable to create shared file %s in openfile\n", rbasciiFilename);
				fflush(stdout);
			}
		}
		else if(ASCII_FLAG == 2) // i.e. io_option = 8, rb M files
		{
			rc = MPI_File_open(MPI_COMM_SELF, rbnmmFilename, MPI_MODE_CREATE | MPI_MODE_RDWR , MPI_INFO_NULL, &mfile);
			if(rc){
				printf("Unable to create unique file %s in openfile in rank %d\n", rbnmmFilename, myrank);
				fflush(stdout);
			}
		}
	}//end of if mySpecies == 1

	if(DEBUG_FLAG) printf("openfile6() done\n");
}

#ifdef UPCASE
void PVTK_NMM(  int *id)
#elif  IBM
void pvtk_nmm(  int *id)
#else
void pvtk_nmm_(  int *id)
#endif
{
	printf("0 getting into pvtk_nmm, *id is %d\n", *id);
	char *pvtkcontent;
	pvtkcontent = (char*)malloc(128 * numGroups );
	memset((void*)pvtkcontent, '\0', 128*numGroups);

	sprintf(pvtkcontent, " <File version=\"pvtk-1.0\"\n dataType=\"vtkUnstructuredGrid\"\n   numberOfPieces=\"%.6d\">\n", numGroups);

	for(i = 0; i < numGroups; i++)
	{
		sprintf(pvtkcontent, "%s   \n<Piece fileName=\"mpi-binary-NMM-p%.6d-t%.5d.vtk\"/>",pvtkcontent, i, *id);
	}
	sprintf(pvtkcontent, "%s\n </File>",pvtkcontent); //FIXME misun 7/5/2012;

	char pvtkfname[128];
	MPI_File mpfile;
	memset((void*)pvtkfname, '\0', 128);
	sprintf(pvtkfname, "./vtk/mpi-binary-NMM-t%.5d.pvtk", *id);
	int rrc = MPI_File_open(MPI_COMM_SELF, pvtkfname, MPI_MODE_CREATE | MPI_MODE_RDWR , MPI_INFO_NULL, &mpfile);
	if(rrc){
		printf("Unable to create unique file %s in openfile\n", rbnmmFilename);
		fflush(stdout);
	}
	MPI_Status write_status;
	MPI_File_write_at(mpfile, 0 , pvtkcontent, strlen(pvtkcontent), MPI_CHAR, &write_status);
	MPI_File_close( & mpfile );
	free(pvtkcontent);
}

#ifdef UPCASE
void PVTK_NM(  int *id)
#elif  IBM
void pvtk_nm(  int *id)
#else
void pvtk_nm_(  int *id)
#endif
{
	printf("0 getting into pvtk_nm, *id is %d\n", *id);
	char *pvtkcontent;
	pvtkcontent = (char*)malloc(128 * numGroups );
	memset((void*)pvtkcontent, '\0', 128*numGroups);

	sprintf(pvtkcontent, " <File version=\"pvtk-1.0\"\n dataType=\"vtkUnstructuredGrid\"\n   numberOfPieces=\"%.6d\">\n", numGroups);

	for(i = 0; i < numGroups; i++)
	{
		sprintf(pvtkcontent, "%s   \n<Piece fileName=\"mpi-binary-NM-p%.6d-t%.5d.vtk\"/>",pvtkcontent, i, *id);
	}
	sprintf(pvtkcontent, "%s\n </File>",pvtkcontent);//FIXME misun 7/5/2012

	char pvtkfname[128];
	MPI_File mpfile;
	memset((void*)pvtkfname, '\0', 128);
	sprintf(pvtkfname, "./vtk/mpi-binary-NM-t%.5d.pvtk", *id);
	int rrc = MPI_File_open(MPI_COMM_SELF, pvtkfname, MPI_MODE_CREATE | MPI_MODE_RDWR , MPI_INFO_NULL, &mpfile);
	if(rrc){
		printf("Unable to create unique file %s in openfile\n", rbnmmFilename);
		fflush(stdout);
	}
	MPI_Status write_status;
	MPI_File_write_at(mpfile, 0 , pvtkcontent, strlen(pvtkcontent), MPI_CHAR, &write_status);
	MPI_File_close( & mpfile );
	free(pvtkcontent);
}

#ifdef UPCASE
void CLOSEFILE6()
#elif  IBM
void closefile6()
#else
void closefile6_()
#endif
{
  if(!THREAD)
    MPI_File_close( & mfile );
//	free_file_struc( file );
	//if(myrank == 0)printf("I/O size is %ld bytes, numGroup is %d\n",mfileCur, numGroups);

  // TODO: if it's threading I/O and this is the last I/O step, maybe try to
  // join threads here so that everything is finished?
  io_time_end = MPI_Wtime();
  io_time = io_time_end - io_time_start; // time from open to close (not including wait lock time)
}

void workersend()
{
	MPI_Barrier(MPI_COMM_WORLD);

	int destrank = -1;
	//compute destrank here
	destrank = groupSize * groupRank ;

	sendBufferCur = 0;

	//if(AUGMENT_FLAG == 1) //do some augment data for validation
	memcpy(&sendBuffer[sendBufferCur], &rankInGroup, sizeof(int));
	sendBufferCur += sizeof(int);
	memcpy(&sendBuffer[sendBufferCur], &mfBufferCur, sizeof(long long));
	sendBufferCur += sizeof(long long);
	memcpy(&sendBuffer[sendBufferCur], mfBuffer, mfBufferCur);
	sendBufferCur += mfBufferCur;


	//added timing info for Isend function itself and calculate perceived speed
	MPI_Request isend_req;
	long long isend_start, isend_end, isend_cycles, isend_size, isend_totalsize, isend_totalcycles, isend_maxcycles;
	double isend_totaltime, isend_maxtime, isend_avgtime;
	isend_size = sendBufferCur;
	MPI_Barrier(localcomm);
	isend_start = rdtsc();
	MPI_Isend(sendBuffer, sendBufferCur, MPI_CHAR, destrank, 1, MPI_COMM_WORLD, &isend_req);
	isend_end = rdtsc();
	MPI_Barrier(localcomm);
	isend_cycles = isend_end - isend_start;
	double isend_time = (double) (isend_cycles/CPU_FREQ);
	MPI_Allreduce(&isend_size, &isend_totalsize,  1, MPI_LONG_LONG_INT, MPI_SUM, localcomm);
	MPI_Allreduce(&isend_cycles, &isend_totalcycles,  1, MPI_LONG_LONG_INT, MPI_SUM, localcomm);
	MPI_Allreduce(&isend_cycles, &isend_maxcycles,  1, MPI_LONG_LONG_INT, MPI_MAX, localcomm);
	//MPI_Allreduce(  &isend_time, &isend_totaltime, 1, MPI_DOUBLE, MPI_SUM, localcomm);
	//MPI_Allreduce(  &isend_time, &isend_maxtime, 1, MPI_DOUBLE, MPI_MAX, localcomm);
	isend_avgtime = isend_totaltime/localsize; //FIXME misun 7/5/2012

	long long isend_avgcycles = isend_totalcycles/localsize;
	if(DEBUG_FLAG && localrank == 0) printf("isend total size is %lld bytes, isend avgtime is %lld cycles, isend maxtime is %lld cycles\n", isend_totalsize, isend_avgcycles, isend_maxcycles);
	if(DEBUG_FLAG) printf("sent size = %d, from rank %d to rank %d\n", sendBufferCur, myrank, destrank);

	MPI_Barrier(MPI_COMM_WORLD);
	mfBufferCur = 0, sendBufferCur = 0;

	if(DEBUG_FLAG) printf("workersend() done - myrank = %d\n", myrank);
}

void flushCurrentBuf8()
{
//#ifdef HYBRID_PTHREAD
  if(THREAD) {
	// copy data to thread writer buffer
	if(DEBUG_FLAG) {
		printf("in flushCurrentBuf8(), going to do memcpy(),wBufCur = %lld, rank = %d\n",writerBufferCur, myrank);
	}
	memcpy(&(file->pwriterBuffer[file->llwriterBufferCur]), writerBuffer, writerBufferCur);
	file->llwriterBufferCur += writerBufferCur;
	writerBufferCur = 0;
	if(DEBUG_FLAG)printf("flushCurrentBuf8() done, rank = %d\n", myrank);
  }
//#else
  else {
	MPI_Status write_status;
	MPI_File_write_at(mfile,
                    mfileCur,
                    writerBuffer,
                    writerBufferCur,
                    MPI_CHAR,
                    &write_status);
	mfileCur += writerBufferCur;
	writerBufferCur = 0;
  }
//#endif
}

void throwToDisk()
{
	long long my_data_offset = 0;

	if(ASCII_FLAG == 0 || ASCII_FLAG == 1) {
		writerBufferCur  = 0;
		for( i = 0; i < groupSize; i++) {
			memcpy(&writerBuffer[writerBufferCur], recvBuffers[i], iSize[i]);
			writerBufferCur += iSize[i];
		}

		MPI_Status write_status;
		MPI_Scan(&writerBufferCur, &my_data_offset, 1, MPI_LONG_LONG_INT, MPI_SUM, localcomm);

		MPI_Allreduce(&writerBufferCur, &fieldSizeSum,  1, MPI_LONG_LONG_INT, MPI_SUM, localcomm);
		my_data_offset += mfileCur;
		my_data_offset -= writerBufferCur;
		MPI_File_write_at_all_begin(mfile, my_data_offset, writerBuffer, writerBufferCur, MPI_CHAR);
		MPI_File_write_at_all_end(mfile, writerBuffer, &write_status);

		mfileCur += fieldSizeSum;

		writerBufferCur = 0;

		if(DEBUG_FLAG)
			printf("throwToDisk(): written size is %lld, file size is %lld\n", fieldSizeSum, mfileCur);
	} // end of ascii_flag = 0 or 1
	else if(ASCII_FLAG == 2) {
		long long thisFieldSize = 0;
		for( i = 0; i < groupSize; i++)
			thisFieldSize += iSize[i];

		if(writerBufferCur + thisFieldSize > writerBufferSize) {
			if(DEBUG_FLAG)printf("buffer not enough, flush some data first, rank = %d\n", myrank);
			flushCurrentBuf8();
		}

		/*memcpy recvBuffers to writerBuffers*/
		for( i = 0; i < groupSize; i++) {
			memcpy(&writerBuffer[writerBufferCur], recvBuffers[i], iSize[i]);
			writerBufferCur += iSize[i];
		}

		/*if I'm the last field, flush to disk and I'm done for the file*/
		if(ifield == numFields) {
			flushCurrentBuf8();
			// run i/o thread here, if HYBRID_PTHREAD not define, then it does not
			// do anything
			if(DEBUG_FLAG)printf("going into run_io_thread(), rank = %d\n", myrank);
			run_io_thread(file);
		}
	} // end of ascii_flag = 2
} // end of throwToDisk

void writerreceive()
{
	MPI_Barrier(MPI_COMM_WORLD);

	/*iSize[i] contains field size (in bytes) of proc with rankIngroup = i*/
	iSize[0] = mfBufferCur;
	memcpy(recvBuffers[0], mfBuffer, mfBufferCur);
	mfBufferCur = 0;

	MPI_Status recv_sta;
	int intrank = -1;
	long long intsize = -1;

	int irecvmsgBufferCur = 0;
	for( i = 1; i < groupSize; i++)
	{
		MPI_Recv(recvmsgBuffer, fieldSizeLimit, MPI_CHAR, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &recv_sta);
		irecvmsgBufferCur = 0;
		memcpy(&intrank, &recvmsgBuffer[irecvmsgBufferCur], sizeof(int));
		irecvmsgBufferCur += sizeof(int);
		memcpy(&intsize, &recvmsgBuffer[irecvmsgBufferCur], sizeof(long long));
		irecvmsgBufferCur += sizeof(long long);
		iSize[intrank] = intsize;

		memcpy(recvBuffers[intrank], &recvmsgBuffer[irecvmsgBufferCur], intsize);
		if(DEBUG_FLAG)printf("writer %d received size = %lld from rank %d \n",myrank, intsize, intrank);
	}

	ifield ++;


	/*write data to disks*/
	if(DEBUG_FLAG)printf("going into throwToDisk() for ifield = %d, rank = %d\n", ifield, myrank);
	throwToDisk();

	MPI_Barrier(MPI_COMM_WORLD);
}

void writeOutField()
{
	if(ASCII_FLAG == 0 || ASCII_FLAG == 1 || ASCII_FLAG == 2)
	{
		if     (mySpecies == 2)
			workersend();
		else if(mySpecies == 1)
			writerreceive();
	}
	/* 3 means it's NM collctive, diff from other NMx case*/
	else if (ASCII_FLAG == 3)
	{
		long long my_data_offset = 0;
		MPI_Status write_status;
		MPI_Scan(&mfBufferCur, &my_data_offset, 1, MPI_LONG_LONG_INT, MPI_SUM, groupcomm);
		MPI_Allreduce(&mfBufferCur, &fieldSizeSum,  1, MPI_LONG_LONG_INT, MPI_SUM, groupcomm);

		my_data_offset += mfileCur;
		my_data_offset -= mfBufferCur;
		MPI_File_write_at_all_begin(mfile, my_data_offset, mfBuffer, mfBufferCur, MPI_CHAR);
		MPI_File_write_at_all_end(mfile, mfBuffer, &write_status);

		mfileCur += fieldSizeSum;
		mfBufferCur = 0;
	}
}

#ifdef UPCASE
void WRITEHEADER6()
#elif  IBM
void writeheader6()
#else
void writeheader6_()
#endif
{
	char* sHeader = (char*)malloc(1024 * sizeof(char));
	memset((void*)sHeader, '\0', 1024);
	sprintf(sHeader, "# vtk DataFile Version 3.0 \n");
	sprintf(sHeader+strlen(sHeader), "Electromagnetic Field  \n");
	if(ASCII_FLAG == 0 || ASCII_FLAG == 2 || ASCII_FLAG == 3)sprintf(sHeader+strlen(sHeader),  "BINARY \n");
	else if(ASCII_FLAG == 1) sprintf(sHeader+strlen(sHeader),  "ASCII \n");
	sprintf(sHeader+strlen(sHeader), "DATASET UNSTRUCTURED_GRID \n");

	if( ((ASCII_FLAG ==0 || ASCII_FLAG == 1) && myrank == 0)  ||   ((ASCII_FLAG == 2 || ASCII_FLAG == 3) && rankInGroup == 0)) //HARDCODED for now, should be first local_rank
	{
		memcpy(mfBuffer, sHeader, strlen(sHeader));
		mfBufferCur += strlen(sHeader);
	}

	free(sHeader);
}

#ifdef UPCASE
void WRITENODES6(double *xyzCoords, int *numNodes)
#elif  IBM
void writenodes6(double *xyzCoords, int *numNodes)
#else
void writenodes6_(double *xyzCoords, int *numNodes)
#endif
{
	float coord[3];
	int   i, j;
	/*
		 fprintf(fp, "POINTS  %d ", *numNodes );
		 fprintf(fp, " float  \n");
		 */
	int totalNumNodes = 0;
	if(ASCII_FLAG == 0 || ASCII_FLAG == 1)
		MPI_Allreduce(numNodes, &totalNumNodes, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD);
	else if(ASCII_FLAG == 2 || ASCII_FLAG == 3)
		MPI_Allreduce(numNodes, &totalNumNodes, 1, MPI_INTEGER, MPI_SUM, groupcomm);
	if( ((ASCII_FLAG ==0 || ASCII_FLAG == 1) && myrank == 0)  ||   ((ASCII_FLAG == 2 ||ASCII_FLAG == 3) && rankInGroup == 0))
	{
		char* sHeader = (char*) malloc( 1024*sizeof(char));
		memset((void*)sHeader, '\0', 1024);
		sprintf(sHeader, "POINTS  %d ", totalNumNodes );
		sprintf(sHeader+strlen(sHeader),  " float  \n");

		memcpy(&mfBuffer[mfBufferCur], sHeader, strlen(sHeader));
		mfBufferCur += strlen(sHeader);
		free(sHeader);
	}

	for( i = 0; i < *numNodes; i++) {
		coord[0] = (float)xyzCoords[3*i+0];
		coord[1] = (float)xyzCoords[3*i+1];
		coord[2] = (float)xyzCoords[3*i+2];

		if(ASCII_FLAG == 0 || ASCII_FLAG == 2 || ASCII_FLAG == 3)
		{
			swap_float_byte( &coord[0] );
			swap_float_byte( &coord[1] );
			swap_float_byte( &coord[2] );

			memcpy(&mfBuffer[mfBufferCur], coord, sizeof(float)*3);
			mfBufferCur += sizeof(float) *3;
		}
		else if( ASCII_FLAG == 1)
		{
			for(j = 0; j < 3; j++)
			{
				sprintf(&mfBuffer[mfBufferCur], "%18.8E", coord[j]);
				mfBufferCur += FLOAT_DIGITS;
			}
			sprintf(&mfBuffer[mfBufferCur++], "\n");
		}
	}

	writeOutField();

	/*fprintf( fp, " \n");*/
	/*simply adding "\n", following original format*/
	if( ((ASCII_FLAG ==0 || ASCII_FLAG == 1) && myrank == 0)  ||   ((ASCII_FLAG == 2 || ASCII_FLAG == 3) && rankInGroup == 0))
	{
		mfBuffer[mfBufferCur] = '\n';
		mfBufferCur++;
	}


}


#ifdef UPCASE
void WRITE2DCELLS6( int *eConnect, int *numElems, int *numCells, int *numNodes)
#elif  IBM
void write2dcells6( int *eConnect, int *numElems, int *numCells, int *numNodes)
#else
void write2dcells6_( int *eConnect, int *numElems, int *numCells, int *numNodes)
#endif
{
	int conn[5];
	int i, j;
	int elemType=9;

	//cell number would be aggregated here
	////following conn number would add an offset - myeConnOffset
	int totalNumCells = 0;
	int myeConnOffset = 0;
	int totalNumNodes = 0;

	if(ASCII_FLAG == 0 || ASCII_FLAG == 1)
	{
		MPI_Allreduce( numCells, &totalNumCells, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD);

		MPI_Scan( numNodes, &myeConnOffset, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD);
		myeConnOffset -= *numNodes;

		MPI_Allreduce(numNodes, &totalNumNodes, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD);
	}
	else if( ASCII_FLAG == 2 || ASCII_FLAG == 3)
	{
		MPI_Allreduce( numCells, &totalNumCells, 1, MPI_INTEGER, MPI_SUM, groupcomm);

		MPI_Scan( numNodes, &myeConnOffset, 1, MPI_INTEGER, MPI_SUM, groupcomm);
		myeConnOffset -= *numNodes;

		MPI_Allreduce(numNodes, &totalNumNodes, 1, MPI_INTEGER, MPI_SUM, groupcomm);
	}

	if( ((ASCII_FLAG ==0 || ASCII_FLAG == 1) && myrank == 0)  ||   ((ASCII_FLAG == 2 || ASCII_FLAG == 3)&& rankInGroup == 0))
	{
		char* sHeader = (char*) malloc (1024 * sizeof(char));
		memset((void*)sHeader, '\0', 1024);
		sprintf(sHeader, "CELLS %d  %d \n", totalNumCells, 5*(totalNumCells) );

		memcpy(&mfBuffer[mfBufferCur], sHeader, strlen(sHeader));
		mfBufferCur += strlen(sHeader);
		free(sHeader);
	}


	for (i = 0; i < *numCells; i++) {

		conn[0] = 4;
		conn[1] = eConnect[4*i+0] + myeConnOffset;
		conn[2] = eConnect[4*i+1] + myeConnOffset;
		conn[3] = eConnect[4*i+2] + myeConnOffset;
		conn[4] = eConnect[4*i+3] + myeConnOffset;

		if(ASCII_FLAG == 0 || ASCII_FLAG == 2 || ASCII_FLAG == 3)
		{
			for( j = 0; j < 5; j++) swap_int_byte( &conn[j] );

			memcpy(&mfBuffer[mfBufferCur], conn, sizeof(int)*5);
			mfBufferCur += sizeof(int) * 5;
		}
		else if( ASCII_FLAG == 1)
		{
			for( j = 0 ; j < 5; j ++)
			{
				sprintf(&mfBuffer[mfBufferCur], "%10d", conn[j]);
				mfBufferCur += INT_DIGITS;
			}
			sprintf(&mfBuffer[mfBufferCur++], "\n");

		}
	}
	//flush to disk

	writeOutField();

	/*
		 fprintf( fp, "\n");
		 fprintf( fp, "CELL_TYPES %d \n", *numCells);
		 */

	if( ((ASCII_FLAG ==0 || ASCII_FLAG == 1) && myrank == 0)  ||   ((ASCII_FLAG == 2 || ASCII_FLAG == 3) && rankInGroup == 0))
	{
		char* sHeader = (char*) malloc (1024 * sizeof(char));
		memset((void*)sHeader, '\0', 1024);
		sprintf(sHeader, "\nCELL_TYPES %d \n", totalNumCells );

		memcpy(&mfBuffer[mfBufferCur], sHeader, strlen(sHeader));
		mfBufferCur += strlen(sHeader);
		free(sHeader);
	}

	//ascii case does not need swapping
	if(ASCII_FLAG != 1)
		swap_int_byte(&elemType);

	for( i = 0; i < *numCells; i++)
	{
		//fwrite(&elemType,  sizeof(int), 1, fp);
		if(ASCII_FLAG == 0 || ASCII_FLAG == 2 || ASCII_FLAG == 3)
		{
			memcpy(&mfBuffer[mfBufferCur], &elemType, sizeof(int));
			mfBufferCur += sizeof(int);
		}
		else if(ASCII_FLAG == 1)
		{
			sprintf(&mfBuffer[mfBufferCur], "%4d", elemType);
			mfBufferCur += 4;

		}
	}

	writeOutField();
	/*
		 fprintf( fp, "\n");
		 fprintf( fp, "POINT_DATA %d \n", *numNodes);
		 */
	if( ((ASCII_FLAG ==0 || ASCII_FLAG == 1) && myrank == 0)  ||
			((ASCII_FLAG == 2 || ASCII_FLAG == 3)&& rankInGroup == 0))
	{
		char* sHeader = (char*) malloc (1024 * sizeof(char));
		memset((void*)sHeader, '\0', 1024);
		sprintf(sHeader, "\nPOINT_DATA %d \n", totalNumNodes );

		memcpy(&mfBuffer[mfBufferCur], sHeader, strlen(sHeader));
		mfBufferCur += strlen(sHeader);
		free(sHeader);
	}
}

#ifdef UPCASE
void WRITE3DCELLS6( int *eConnect, int *numElems, int *numCells, int *numNodes)
#elif  IBM
void write3dcells6( int *eConnect, int *numElems, int *numCells, int *numNodes)
#else
void write3dcells6_( int *eConnect, int *numElems, int *numCells, int *numNodes)
#endif
{
	int conn[9];
	int conn_new[9];
	int i, j;
	int elemType=12;

	int totalNumCells = 0;
	int myeConnOffset = 0;
	int totalNumNodes = 0;

	if(ASCII_FLAG == 0 || ASCII_FLAG == 1)
	{
		MPI_Allreduce( numCells, &totalNumCells, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD);

		MPI_Scan( numNodes, &myeConnOffset, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD);
		myeConnOffset -= *numNodes;

		MPI_Allreduce(numNodes, &totalNumNodes, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD);
	}
	else if( ASCII_FLAG == 2 || ASCII_FLAG == 3)
	{
		MPI_Allreduce( numCells, &totalNumCells, 1, MPI_INTEGER, MPI_SUM, groupcomm);

		MPI_Scan( numNodes, &myeConnOffset, 1, MPI_INTEGER, MPI_SUM, groupcomm);
		myeConnOffset -= *numNodes;

		MPI_Allreduce(numNodes, &totalNumNodes, 1, MPI_INTEGER, MPI_SUM, groupcomm);
	}

	if( ((ASCII_FLAG ==0 || ASCII_FLAG == 1) && myrank == 0)  ||   ((ASCII_FLAG == 2 || ASCII_FLAG == 3) && rankInGroup == 0))
	{
		char* sHeader = (char*) malloc (1024 * sizeof(char));
		memset((void*)sHeader, '\0', 1024);
		sprintf(sHeader, "CELLS %d  %d \n", totalNumCells, 9*(totalNumCells) );

		memcpy(&mfBuffer[mfBufferCur], sHeader, strlen(sHeader));
		mfBufferCur += strlen(sHeader);
		free(sHeader);
	}

	for (i = 0; i < *numCells; i++) {
		/*
			 conn[0] = 8;
			 conn[1] = eConnect[8*i+0];
			 conn[2] = eConnect[8*i+1];
			 conn[3] = eConnect[8*i+2];
			 conn[4] = eConnect[8*i+3];
			 conn[5] = eConnect[8*i+4];
			 conn[6] = eConnect[8*i+5];
			 conn[7] = eConnect[8*i+6];
			 conn[8] = eConnect[8*i+7];
			 for( j = 0; j < 9; j++) swap_int_byte( &conn[j] );
			 fwrite(conn, sizeof(int), 9, fp);
			 */

		conn_new[0] = 8;
		conn_new[1] = eConnect[8*i+0] + myeConnOffset;
		conn_new[2] = eConnect[8*i+1] + myeConnOffset;
		conn_new[3] = eConnect[8*i+2] + myeConnOffset;
		conn_new[4] = eConnect[8*i+3] + myeConnOffset;
		conn_new[5] = eConnect[8*i+4] + myeConnOffset;
		conn_new[6] = eConnect[8*i+5] + myeConnOffset;
		conn_new[7] = eConnect[8*i+6] + myeConnOffset;
		conn_new[8] = eConnect[8*i+7] + myeConnOffset;

		if(ASCII_FLAG == 0 || ASCII_FLAG == 2 || ASCII_FLAG == 3)
		{
			for( j = 0; j < 9; j++) swap_int_byte( &conn_new[j] );

			memcpy(&mfBuffer[mfBufferCur], conn_new, sizeof(int)*9);
			mfBufferCur += sizeof(int) * 9;
		}
		else if( ASCII_FLAG == 1)
		{
			for( j = 0; j < 9; j ++)
			{
				sprintf(&mfBuffer[mfBufferCur], "%10d", conn_new[j]);
				mfBufferCur += INT_DIGITS ;
			}
			sprintf(&mfBuffer[mfBufferCur++], "\n");

		}
	}
	writeOutField();

	/*
		 fprintf( fp, "\n");
		 fprintf( fp, "CELL_TYPES %d \n", *numCells );
		 */
	if( ((ASCII_FLAG ==0 || ASCII_FLAG == 1) && myrank == 0)  ||   ((ASCII_FLAG == 2 || ASCII_FLAG == 3) && rankInGroup == 0))
	{
		char* sHeader = (char*) malloc (1024 * sizeof(char));
		memset((void*)sHeader, '\0', 1024);
		sprintf(sHeader, "\nCELL_TYPES %d \n", totalNumCells );

		memcpy(&mfBuffer[mfBufferCur], sHeader, strlen(sHeader));
		mfBufferCur += strlen(sHeader);
		free(sHeader);
	}

	// ascii case does not need swapping (ascii is universal format)
	if(ASCII_FLAG != 1)
		swap_int_byte(&elemType);

	for (i = 0; i < *numCells; i++)
	{
		//fwrite(&elemType,  sizeof(int), 1, fp);

		if(ASCII_FLAG == 0 || ASCII_FLAG == 2 || ASCII_FLAG == 3)
		{
			memcpy(&mfBuffer[mfBufferCur], &elemType, sizeof(int));
			mfBufferCur += sizeof(int);
		}
		else if( ASCII_FLAG == 1)
		{
			sprintf(&mfBuffer[mfBufferCur], "%4d", elemType);
			mfBufferCur += 4 ;

		}

	}
	writeOutField();

	/*
		 fprintf( fp, "\n");
		 fprintf( fp, "POINT_DATA %d \n", *numNodes );
		 */
	if( ((ASCII_FLAG ==0 || ASCII_FLAG == 1) && myrank == 0)  ||   ((ASCII_FLAG == 2 || ASCII_FLAG == 3) && rankInGroup == 0))
	{
		char* sHeader = (char*) malloc (1024 * sizeof(char));
		memset((void*)sHeader, '\0', 1024);
		sprintf(sHeader, "\nPOINT_DATA %d \n", totalNumNodes );

		memcpy(&mfBuffer[mfBufferCur], sHeader, strlen(sHeader));
		mfBufferCur += strlen(sHeader);
		free(sHeader);
	}

}

#ifdef UPCASE
void WRITEFIELD6(int *fldid, double *vals, int *numNodes)
#elif  IBM
void writefield6(int *fldid, double *vals, int *numNodes)
#else
void writefield6_(int *fldid, double *vals, int *numNodes)
#endif
{
	float fldval[3];
	int   i, j  ;
	char  fldname[100];

	getfieldname_(*fldid, fldname);
	/*
		 fprintf( fp, "VECTORS %s ", fldname);
		 fprintf( fp, " float \n");
		 */

	if( ((ASCII_FLAG ==0 || ASCII_FLAG == 1) && myrank == 0)  ||   ((ASCII_FLAG == 2 || ASCII_FLAG == 3) && rankInGroup == 0))
	{
		char* sHeader = (char*) malloc (1024 * sizeof(char));
		memset((void*)sHeader, '\0', 1024);
		sprintf(sHeader, "VECTORS %s  float \n", fldname);

		memcpy(&mfBuffer[mfBufferCur], sHeader, strlen(sHeader));
		mfBufferCur += strlen(sHeader);
		free(sHeader);
	}

	for (i = 0; i < *numNodes; i++) {

		fldval[0] = (float)vals[3*i+0];
		fldval[1] = (float)vals[3*i+1];
		fldval[2] = (float)vals[3*i+2];
		//fwrite(fldval, sizeof(float), 3, fp);

		if(ASCII_FLAG == 0 || ASCII_FLAG == 2 || ASCII_FLAG == 3)
		{
			swap_float_byte( &fldval[0]);
			swap_float_byte( &fldval[1]);
			swap_float_byte( &fldval[2]);

			memcpy( &mfBuffer[mfBufferCur], fldval, sizeof(float)*3);
			mfBufferCur += sizeof(float) *3;
		}
		else if( ASCII_FLAG == 1)
		{
			for( j = 0; j < 3; j++)
			{
				sprintf(&mfBuffer[mfBufferCur], "%18.8E", fldval[j]);
				mfBufferCur += FLOAT_DIGITS;
			}
			sprintf(&mfBuffer[mfBufferCur++], "\n");

		}
	}
	writeOutField();

	//   fprintf(fp, " \n");
	//add this return symbol into mpifile...
	if( ((ASCII_FLAG ==0 || ASCII_FLAG == 1) && myrank == 0) ||
			((ASCII_FLAG == 2 || ASCII_FLAG == 3) && rankInGroup == 0))
	{
		char* sHeader = (char*) malloc (1024 * sizeof(char));
		memset((void*)sHeader, '\0', 1024);
		sprintf(sHeader, " \n");

		memcpy(&mfBuffer[mfBufferCur], sHeader, strlen(sHeader));
		mfBufferCur += strlen(sHeader);
		free(sHeader);
	}
}

/**
 * This function allocates and initiates necessary components of file struct, e.g. thread
 * This function is called once per program life time
 *
 */
void init_file_struc() {
//#ifdef HYBRID_PTHREAD
  if(THREAD) {
	if(myrank == 0 && DEBUG_FLAG)printf("Hybrid_pthread defined, using thread rbIO. Now initializing...\n");
	if (file == NULL) {
		if(DEBUG_FLAG) printf("init file_struc for first time\n");
		file = (file_t*) malloc (sizeof(file_t));
		// TODO: accurately estimate buffer size needed
    // Can't re-use the buffer alloc'ed from writerBufferSize coz' this buffer
    // is going to live even after the computation is finished
		file->pwriterBuffer = (char*) malloc( sizeof(char) * WRITERBUFFERSIZE);
		assert(file->pwriterBuffer != NULL);
		memset(file->pwriterBuffer, '\0', WRITERBUFFERSIZE);
		file->pmfile = &mfile;

		pthread_mutex_init(&file->mutex, NULL);
		pthread_cond_init(&file->cond, NULL);
	}
//#endif
  }
}

/**
 * This function resets the file struc, i.e. re-claim the lock and reset the
 * buffer position
 * This function is called once per checkpoint step
 */
void reset_file_struc(){
//#ifdef HYBRID_PTHREAD
  if(THREAD) {
	// if trylock() returns busy, will print a message and block wait
	// otherwise we would have got the lock
	if( pthread_mutex_trylock(&file->mutex) == EBUSY) {
		if(myrank == 0) printf("WARNING: there is an I/O thread grabing the lock.. waiting..\n");
    if( strstr(mach_name, "Intrepid") != NULL && io_step == 1)
      printf("This run is on Intrepid, did you remember to submit job in CO/SMP mode?\n");
		// blocking wait
		pthread_mutex_lock(&file->mutex);
	}
  else {
    if(DEBUG_FLAG) printf("pthread_mutex_trylock succeeded\n");
  }
	file->llwriterBufferCur = 0;
  }
//#endif
}

/**
 * This function frees file struct
 *
 * @param	file	pointer to the file struct
 */
void free_file_struc(file_t* file) {
//#ifdef HYBRID_PTHREAD
  if(THREAD) {
	free(file->pwriterBuffer);

	pthread_mutex_destroy(&file->mutex);
	pthread_cond_destroy(&file->cond);

	free(file);
  }
//#endif
}

/** This function simply writes whole buffer into file using COMM_SELF
 *
 * @param file  pointer to the file struct
 */
void* write_file_buffer(void* arg) {
	file_t* myfile = (file_t*) arg;

	if(myrank == 0 && DEBUG_FLAG) printf("write_file_buffer() - to write %lld bytes\n", myfile->llwriterBufferCur);
	double startio, endio;
	startio = MPI_Wtime();
  // using diff options to see the effect of I/O thread v.s. computation
  int option = 1;
  if(option == 1) {
	MPI_Status write_status;
	MPI_File_write_at(*(myfile->pmfile), 0, myfile->pwriterBuffer,
									  myfile->llwriterBufferCur, MPI_CHAR, &write_status);
  }
  else if (option == 2) {
    int fd = open(rbnmmFilename, O_WRONLY | O_CREAT);
    if (fd == -1) {
      printf("Error: File %s can't be opened (in POSIX)\n", rbnmmFilename);
      exit(1);
    }
    write(fd, myfile->pwriterBuffer, (size_t) myfile->llwriterBufferCur);
    close(fd);
    if(myrank == 0) printf("Using POSIX write\n");
  }
  else if (option == 3) {
    int sleep_sec = 20;
    sleep(sleep_sec);
    if(myrank == 0) printf("Using sleep %d sec\n", sleep_sec);
  }
	endio = MPI_Wtime();
  file_io_time = endio - startio;

	MPI_File_close( myfile->pmfile );

	if( myrank == 0 && DEBUG_FLAG)
    printf("\nINFO:io thread finished writing one file..took %f sec, file_io_time = %lf sec- rank = %d\n", endio - startio, file_io_time, myrank);
  pthread_mutex_unlock(&file->mutex);
	pthread_exit(NULL);
}

/**
 * This function creates io thread, write the data, and exit the pthread
 * This function is called once for every file write
 *
 * @param file  pointer to the file struct
 */
void run_io_thread(file_t* file){
//#ifdef HYBRID_PTHREAD
  if(THREAD) {
	int rc;
	//rc = pthread_create(&file->pthread, NULL, void *(*write_file_buffer) (void*), (void*)file); //FIXME misun 7/5/1202
	rc = pthread_create(&file->pthread, NULL, (void*)write_file_buffer, (void*)file);
	if(rc) {
		printf("ERROR: pthread_create() failed, code: %d\n", rc);
    if( strstr(mach_name, "Intrepid") != NULL )
      printf("Note: this run is on Intrepid, did you remember to submit job in CO/SMP mode?\n");
		exit(-1);
	}
	// one can't call pthread_exit here because  it's main and will hang all
	// following thread
	//pthread_exit(NULL);
  }
//#endif
}
