#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <stdint.h>

#include <mpi.h>
#include "rdtsc.h"
#include "vtkcommon.h"

#define WRITERBUFFERSIZE (50 * ONE_MILLION)

extern MPI_File mfile;
int64_t testint64;
int fieldSizeLimit ;

char rbFilename[100];
char rbasciiFilename[100];
//char mfBuffer[ 4* ONE_MILLION];
//long long mfBufferCur = 0;
//long long mfileCur = 0;
//long long fieldSizeSum = 0;

char *sendBuffer;
int sendBufferCur = 0;

/*global rank and size*/
int myrank, mysize; 
//rank and size in localcomm
int localrank, localsize; 
//worker comm and write comm 
MPI_Comm localcomm, groupcomm; 
int numGroups, groupSize, rankInGroup, groupRank,  mySpecies, numFields, ifield;

char* writerBuffer;
char** recvBuffers;
char* recvmsgBuffer;
int* iSize;
long long writerBufferCur = 0;
int writerBufferSize ;
int recvmsgBufferCur = 0;

int first_init = 0;
int AUGMENT_FLAG = 0;
int ASCII_FLAG = 0; 

int INT_DIGITS = 10;
int FLOAT_DIGITS = 18;
int LONG_LONG_DIGITS = 18;

#ifdef UPCASE
void SET_ASCII_TRUE()
#elif  IBM
void set_ascii_true()
#else
void set_ascii_true_()
#endif
{
	ASCII_FLAG = 1;
	if(DEBUG_FLAG)printf("setting ascii flag to true");
}

#ifdef UPCASE
void SET_ASCII_NMM()
#elif  IBM
void set_ascii_nmm()
#else
void set_ascii_nmm_()
#endif
{
        ASCII_FLAG = 2;
        if(DEBUG_FLAG)printf("setting ascii flag to NMM");
}

#ifdef UPCASE
void SET_ASCII_NM()
#elif  IBM
void set_ascii_nm()
#else
void set_ascii_nm_()
#endif
{
        ASCII_FLAG = 3;
        if(DEBUG_FLAG)printf("setting ascii flag to NM");
}

#ifdef UPCASE
void FREE_RBIO_BUFFER()
#elif  IBM
void free_rbio_buffer()
#else
void free_rbio_buffer_()
#endif
{
	free(mfBuffer);
	free(sendBuffer);
	for( int i = 0; i < groupSize; i++)
                        free (recvBuffers[i]) ;
	free(recvBuffers);
	free(iSize);
	free(recvmsgBuffer);
	free(writerBuffer);
}


void smartCheckGroupSize(int *numgroups)
{
	int IDEAL_SIZE = 64;
	int SIZE_UPPER_BOUND = 1024;
	if( (mysize/(*numgroups)) > SIZE_UPPER_BOUND)
	{
		*numgroups = mysize/IDEAL_SIZE;
		if(myrank == 0)printf("changed numGroups to %d\n", *numgroups);
	}
}

#ifdef UPCASE
void INITRBIO(int *numgroups, int* numfields, int* maxnumnodes)
#elif  IBM
void initrbio(int *numgroups, int* numfields, int* maxnumnodes)
#else
void initrbio_(int *numgroups, int* numfields, int* maxnumnodes)
#endif
{   if (first_init == 0)
	{
	first_init = 1; //only init for first time

	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &mysize);

	smartCheckGroupSize(numgroups);

	numGroups = *numgroups;
	//printf("numgroups is %d*****\n", *numgroups);
	numFields = *numfields + 3; //+3 due to node coord, cell numbering and celltype 	

	groupSize = mysize / numGroups;
	groupRank = myrank / groupSize;
	rankInGroup = myrank % groupSize;

	/*
 	*upper bound of single field size(include nodes info and 3d cells)
	*for ascii case, float has 18 digits * 3, int has 7 digits * 10, so take int size
	*/
	if(ASCII_FLAG == 0 || ASCII_FLAG == 2 || ASCII_FLAG == 3)
	fieldSizeLimit = 10 * sizeof(int) * (*maxnumnodes) + 1024;
	else if(ASCII_FLAG == 1)
	fieldSizeLimit = 10 * (INT_DIGITS + 1) * (*maxnumnodes) + 1024;

	if(ASCII_FLAG == 0 || ASCII_FLAG == 1 || ASCII_FLAG == 3)
	writerBufferSize = groupSize * fieldSizeLimit;
	else if(ASCII_FLAG == 2)
	writerBufferSize = WRITERBUFFERSIZE;	
	
	mfBuffer = (char*) malloc( sizeof(char) * fieldSizeLimit);
	sendBuffer = (char*) malloc( sizeof(char) * fieldSizeLimit);

	//writer species is 1, worker is 2
	if(rankInGroup == 0) 
		{
		mySpecies = 1; 
		recvBuffers = (char**)malloc(sizeof(char*) * groupSize);
		for( int i = 0; i < groupSize; i++)
			recvBuffers[i] = (char*) malloc(sizeof(char) * fieldSizeLimit); 
		iSize = (int*) malloc(sizeof(int) * groupSize);
		recvmsgBuffer = (char*) malloc( sizeof(char) * fieldSizeLimit);
		writerBuffer = (char*) malloc(sizeof(char) * writerBufferSize);
		if(writerBuffer == NULL)
		{
			printf("not enough memory for writerBuffer - quitting!!\n");
			exit(1);
		}
		}	
	else mySpecies = 2;
	
	MPI_Comm_split(MPI_COMM_WORLD, mySpecies, myrank, &localcomm);
	MPI_Comm_split(MPI_COMM_WORLD, groupRank, myrank, &groupcomm);
	MPI_Comm_rank(localcomm, &localrank);
	MPI_Comm_size(localcomm, &localsize);
	if(DEBUG_FLAG)printf("myrank is %d, rankInGroup = %d, localrank is %d, numGroups is %d, fieldSizeLimit is %d, maxnumnodes is %d\n", myrank, rankInGroup, localrank, numGroups, fieldSizeLimit, *maxnumnodes);
	}
	
	ifield = 0;
	mfBufferCur = 0;
   	mfileCur = 0; //reset file position for new file
}

#ifdef UPCASE
void OPENFILE6(  int *id, int *nid)
#elif  IBM
void openfile6(  int *id, int *nid)
#else
void openfile6_(  int *id, int *nid)
#endif
{
	if(IOTIMER_FLAG)	
	start_time = rdtsc();

	getfilename_(id, nid);
	
	if(ASCII_FLAG == 3)
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
	if(ASCII_FLAG == 0)
	{
        rc = MPI_File_open(localcomm, rbFilename, MPI_MODE_CREATE | MPI_MODE_RDWR , MPI_INFO_NULL, &mfile);
        if(rc){
            printf("Unable to create shared file %s in openfile\n", rbFilename);
            fflush(stdout);
          }
	}
	else if(ASCII_FLAG == 1)
	{
	rc = MPI_File_open(localcomm, rbasciiFilename, MPI_MODE_CREATE | MPI_MODE_RDWR , MPI_INFO_NULL, &mfile);
        if(rc){
            printf("Unable to create shared file %s in openfile\n", rbasciiFilename);
            fflush(stdout);
          }
	}
	else if(ASCII_FLAG == 2)
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

	for(int i = 0; i < numGroups; i++)
	{
		sprintf(pvtkcontent, "%s   \n<Piece fileName=\"mpi-binary-NMM-p%.6d-t%.5d.vtk\"/>",pvtkcontent, i, *id);
	} 
	sprintf(pvtkcontent, "%s\n </File>");
	
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

        for(int i = 0; i < numGroups; i++)
        {
                sprintf(pvtkcontent, "%s   \n<Piece fileName=\"mpi-binary-NM-p%.6d-t%.5d.vtk\"/>",pvtkcontent, i, *id);
        }
        sprintf(pvtkcontent, "%s\n </File>");

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
   MPI_File_close( & mfile );
	if(myrank == 0)printf("I/O size is %ld bytes, numGroup is %d\n",mfileCur, numGroups); 
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
	double isend_time = (double) (isend_cycles/BGP_FREQ); 
	MPI_Allreduce(&isend_size, &isend_totalsize,  1, MPI_LONG_LONG_INT, MPI_SUM, localcomm);	
MPI_Allreduce(&isend_cycles, &isend_totalcycles,  1, MPI_LONG_LONG_INT, MPI_SUM, localcomm);
MPI_Allreduce(&isend_cycles, &isend_maxcycles,  1, MPI_LONG_LONG_INT, MPI_MAX, localcomm);
	//MPI_Allreduce(  &isend_time, &isend_totaltime, 1, MPI_DOUBLE, MPI_SUM, localcomm);
	//MPI_Allreduce(  &isend_time, &isend_maxtime, 1, MPI_DOUBLE, MPI_MAX, localcomm);
	isend_avgtime = isend_totaltime/localsize;
	long long isend_avgcycles = isend_totalcycles/localsize;
	if(localrank == 0)printf("isend total size is %lld bytes, isend avgtime is %lld cycles, isend maxtime is %lld cycles\n", isend_totalsize, isend_avgcycles, isend_maxcycles);	
	if(DEBUG_FLAG)printf("sent size = %d, from rank %d to rank %d\n", sendBufferCur, myrank, destrank);
	
	MPI_Barrier(MPI_COMM_WORLD);
	mfBufferCur = 0, sendBufferCur = 0;
}

void flushCurrentBuf8()
{
	MPI_Status write_status;
	MPI_File_write_at(mfile, mfileCur, writerBuffer, writerBufferCur, MPI_CHAR, &write_status);
	mfileCur += writerBufferCur;
	writerBufferCur = 0;
}

void throwToDisk()
{
        long long my_data_offset = 0;
	
	if(ASCII_FLAG == 0 || ASCII_FLAG == 1)
	{

		writerBufferCur  = 0;
        	for( int i = 0; i < groupSize; i++)
        	{
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
	}
	else if(ASCII_FLAG == 2)
	{
		long long thisFieldSize = 0;
		for( int i = 0; i < groupSize; i++)
			thisFieldSize += iSize[i];
		
		if(writerBufferCur + thisFieldSize > writerBufferSize)
		{
			flushCurrentBuf8();
		}
		
		/*memcpy recvBuffers to writerBuffers*/
		for( int i = 0; i < groupSize; i++)
                {
                        memcpy(&writerBuffer[writerBufferCur], recvBuffers[i], iSize[i]);
                        writerBufferCur += iSize[i];
                }
		
		
		/*if I'm the last field, flush to disk and I'm done for the file*/
		if(ifield == numFields)
		{
			flushCurrentBuf8();
		}
	}	

}

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
	for( int i = 1; i < groupSize; i++)
	{
		MPI_Recv(recvmsgBuffer, fieldSizeLimit, MPI_CHAR, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &recv_sta);
		irecvmsgBufferCur = 0;
		memcpy(&intrank, &recvmsgBuffer[irecvmsgBufferCur], sizeof(int));
		irecvmsgBufferCur += sizeof(int);
		memcpy(&intsize, &recvmsgBuffer[irecvmsgBufferCur], sizeof(long long));
		irecvmsgBufferCur += sizeof(long long);
		iSize[intrank] = intsize;

		memcpy(recvBuffers[intrank], &recvmsgBuffer[irecvmsgBufferCur], intsize);
		if(DEBUG_FLAG)printf("writer %d received size = %lld from rank %d ",myrank, intsize, intrank);
	}	

	ifield ++;


	/*write data to disks*/
	throwToDisk();
		
	MPI_Barrier(MPI_COMM_WORLD);
	
}

void writeOutField()
{
	if(ASCII_FLAG == 0 || ASCII_FLAG == 1 || ASCII_FLAG == 2)
	{
        	if(mySpecies == 2)
                	workersend();
        	else if(mySpecies == 1)
                	writerreceive();
	}
	/* 3 means it's NM collctive, diff from other NMx case*/
	else if(ASCII_FLAG == 3)
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
   memset((void*)sHeader, "\0", 1024);
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
                memset((void*)sHeader, "\0", 1024);
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
		for(int j = 0; j < 3; j++)
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
		for( int j = 0 ; j < 5; j ++)
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
	if( ((ASCII_FLAG ==0 || ASCII_FLAG == 1) && myrank == 0)  ||   ((ASCII_FLAG == 2 || ASCII_FLAG == 3)&& rankInGroup == 0)) 
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
		for( int j = 0; j < 9; j ++)
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
		for( int j = 0; j < 3; j++)
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
        if( ((ASCII_FLAG ==0 || ASCII_FLAG == 1) && myrank == 0)  ||   ((ASCII_FLAG == 2 || ASCII_FLAG == 3) && rankInGroup == 0)) 
        {
        char* sHeader = (char*) malloc (1024 * sizeof(char));
        memset((void*)sHeader, '\0', 1024);
        sprintf(sHeader, " \n");

        memcpy(&mfBuffer[mfBufferCur], sHeader, strlen(sHeader));
        mfBufferCur += strlen(sHeader);
        free(sHeader);
        }
}
