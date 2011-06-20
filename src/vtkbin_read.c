#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include <mpi.h>
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

#define ONE_MILLION 1048576
//FILE *fp = NULL;

MPI_File mfile;

//char filename[100];
//char mFilename[100];

//char* mfBuffer ;
//long long mfileCur = 0, mfBufferCur = 0;
//long long fieldSizeSum = 0;

//int myrank;

/************************************************
 *  MPI-IO format (param(102)=4, syncIO ASCII)
 ************************************************
 */

/**
 * This function gets proper filename and opens it
 *
 * @param id	checkpoint time step
 * @param nid	processor id
 *
 */

#ifdef UPCASE
void OPENFILE_READ4(  int *id, int *nid)
#elif  IBM
void openfile_read4(  int *id, int *nid)
#else
void openfile_read4_(  int *id, int *nid)
#endif
{
	// to figure out which restart step (thus path) to read
	// should get this value from .rea input file
	// hack it for now to work with local case
	// TODO -- local vtk will be deleted when re-run?
	getfilename_(id,nid, 4);

	/* parallel here*/

	int rc;
	rc = MPI_File_open(MPI_COMM_WORLD, mFilename, MPI_MODE_RDONLY , MPI_INFO_NULL, &mfile);
	if(rc){
		printf("Unable to open shared file (%s) in openfile() to read\n", mFilename);
		fflush(stdout);
	}

	mfBuffer = (char*) malloc( sizeof( char) * 4 * ONE_MILLION);
	mfBufferCur = 0;
}

#ifdef UPCASE
void CLOSEFILE_READ4()
#elif  IBM
void closefile_read4()
#else
void closefile_read4_()
#endif
{
   MPI_File_close( & mfile );
}


#ifdef UPCASE
void READHEADER4(int* irsttemp, int* idump)
#elif  IBM
void readheader4(int* irsttemp, int* idump)
#else
void readheader4_(int* irsttemp, int* idump)
#endif
{
   int i ;
   float xyz[3];

	 /*put header into sHeader string */
   char sHeader[1024];
   memset((void*)sHeader, '\0', 1024);

   mfileCur = 0; //reset file position for new file
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank );
   if(myrank == 0) // need only one proc to read header
   {
		 MPI_Status status;
		 //TODO: how to check status and return value for reliability?
		 // return value = 0 indicates success
		 // sizeof(str) -1 here because sizeof() returns size of string including
		 // extra null, but null is not written into file
		 int ierr = MPI_File_read_at(mfile, mfileCur, sHeader, 1024, MPI_CHAR, &status);
		 if (ierr != 0) {
			 printf("error: MPI_File_read_at failed in readheader4()\n");
			 exit(4);
		 }

		 int headerLen = 0;
		 char* version = strtok(sHeader, "\n");
		 headerLen += strlen(version) + 1; // +1 is to add newline
		 char* restart = strtok(NULL, "\n");
		 headerLen += strlen(restart) + 1;
		 char* ascii = strtok(NULL, "\n");
		 headerLen += strlen(ascii) + 1;
		 char* grid = strtok(NULL, "\n");
		 headerLen += strlen(grid) + 1;
		 //printf("ver: %s, restart: %s, ascii: %s, grid: %s\n", version, restart, ascii, grid);

		 // they are istep, idumpno, time, dt
		 // TODO: error checking here??
		 char* strstep = strtok(restart, " ");
		 *irsttemp = atoi(strstep);
		 char* strdump = strtok(NULL, " ");
		 *idump = atoi(strdump);

		 printf("istep %d idump %d\n", *irsttemp, *idump);
		 mfileCur = headerLen;
   }
	 //if(myrank == 0) printf("readheader4() done\n");
   MPI_Bcast(irsttemp, 1, MPI_INT, 0, MPI_COMM_WORLD); 
   MPI_Bcast(idump, 1, MPI_INT, 0, MPI_COMM_WORLD);  

}

/**
 * This function parse the data header that has only two parts, e.g. cell
 * type, points, etc.
 *
 */
void parseHeader4short(int* len, char* field, char* number) {
	if( myrank == 0) {
		char sHeader[1024];
		memset((void*)sHeader, '\0', 1024);
		int ierr = MPI_File_read_at(mfile, mfileCur, (void*)sHeader, 1024, MPI_CHAR, MPI_STATUSES_IGNORE);
		if (ierr != 0) {
			printf("error: MPI_File_read_at failed in readheader4()\n");
			exit(4);
		}
		//printf("raw string from file read (field, number) is:%s\n", sHeader);

		// now parse the string we have to three fields
		char* pch;
		// find '\n', before newline it's the real header
		pch = strtok(sHeader, "\n");
		*len = strlen(pch);
		*len = *len + 1; // +1 because strtok trimed '\n', which counts as part of header
		char* f;
		char* n;
		f = strtok(pch, " ");
		n = strtok(NULL, " ");
		memcpy(field, f, strlen(f));
		memcpy(number, n, strlen(n));
	}
}

/**
 * This function parse the data header that has three parts, e.g. nodes
 *
 */
void parseHeader4(int* len, char* field, char* number, char* datatype) {
	if( myrank == 0) {
		char sHeader[1024];
		memset((void*)sHeader, '\0', 1024);
		int ierr = MPI_File_read_at(mfile, mfileCur, (void*)sHeader, 1024, MPI_CHAR, MPI_STATUSES_IGNORE);
		if (ierr != 0) {
			printf("error: MPI_File_read_at failed in readheader4()\n");
			exit(4);
		}
		//printf("raw string from file read (field, number, type?) is:%s\n", sHeader);

		// now parse the string we have to three fields
		char* pch;
		// find '\n', before newline it's the real header
		pch = strtok(sHeader, "\n");
		*len = strlen(pch);
		*len = *len + 1; // +1 because strtok trimed '\n', which counts as part of header
		char* f;
		char* n;
		char* d;
		f = strtok(pch, " ");
		n = strtok(NULL, " ");
		d = strtok(NULL, " ");
		memcpy(field, f, strlen(f));
		memcpy(number, n, strlen(n));
		memcpy(datatype, d, strlen(d));
	}
}

#ifdef UPCASE
void READNODES4(double *xyzCoords, int *numNodes)
#elif  IBM
void readnodes4(double *xyzCoords, int *numNodes)
#else
void readnodes4_(double *xyzCoords, int *numNodes)
#endif
{
   float coord[3];
   int   i, j;

	//xyzCoords itself doesn't change (absolute value)
	 long long parsedTotalNumNodes = 0;
	 long long totalNumNodes = 0;

	 //	MPI_Allreduce(numNodes, &totalNumNodes, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD);
	 // rank 0 will read the header and try to figure out number of nodes
	 char field[1024], number[1024], datatype[1024];
	 int headerLen = 0;
	 memset((void*)field, '\0', 1024);
	 memset((void*)number, '\0', 1024);
	 memset((void*)datatype, '\0', 1024);
	 parseHeader4(&headerLen, field, number, datatype);
	 if( myrank == 0)
	 {
		 printf("field is %s, number is %s, datatype is %s, headerLen is %d\n", field, number, datatype, headerLen);
		 // skip header part
		 mfileCur += headerLen;
		 //TODO: what if it's beyond int range?
		 parsedTotalNumNodes = (int) atoi(number);
	 }
	 // rank 0 bcast file pointer position to everyone else
	 MPI_Bcast(&mfileCur, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);

	long long llNumNodes = *numNodes;
	long long myOffset;
	MPI_Scan(&llNumNodes, &myOffset, 1, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);

	// double-check if the parsed number is right
	MPI_Bcast(&parsedTotalNumNodes, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
	if(myrank == mysize && parsedTotalNumNodes != myOffset + (*numNodes)) {
		printf("error: the parsed num nodes is not consistent with solver's \n");
		exit(1);
	}

	float* readxyzCoords= (float*) malloc((*numNodes)*3 *sizeof(float));
	// convert number of elements to number of float numbers
	myOffset -= (*numNodes);
	myOffset *= sizeof(float) *3;

	int ierr = MPI_File_read_at_all( mfile, mfileCur + myOffset, (void*)readxyzCoords, (*numNodes)* sizeof(float) *3, MPI_CHAR, MPI_STATUSES_IGNORE);
	if( ierr != 0) {
		printf("error: MPI_File_read_at_all failed in readnodes4()\n");
		exit(4);
	}
	//update file pointer position
	mfileCur += parsedTotalNumNodes * sizeof(float) * 3;
	// +1 because there is an extra '\n' at end of data block
	mfileCur ++;

   for( i = 0; i < *numNodes; i++) {
		 swap_float_byte( &readxyzCoords[3*i+0]);
		 swap_float_byte( &readxyzCoords[3*i+1]);
		 swap_float_byte( &readxyzCoords[3*i+2]);
		 //printf("read data: %f %f %f\n", readxyzCoords[3*i+0], readxyzCoords[3*i+1], readxyzCoords[3*i+2]);

		 // now pass back to solver
		 xyzCoords[3*i+0] = readxyzCoords[3*i+0];
		 xyzCoords[3*i+1] = readxyzCoords[3*i+1];
		 xyzCoords[3*i+2] = readxyzCoords[3*i+2];
	 }
}

#ifdef UPCASE
void READ2DCELLS4( int *eConnect, int *numElems, int *numCells, int *numNodes)
#elif  IBM
void read2dcells4( int *eConnect, int *numElems, int *numCells, int *numNodes)
#else
void read2dcells4_( int *eConnect, int *numElems, int *numCells, int *numNodes)
#endif
{
   int conn[5];
   int conn_new[5];
   int i, j;
   int elemType=9;

//	 long long totalNumCells = 0;
	 long long parsedTotalNumCells = 0;

	 char field[1024], number[1024], datatype[1024];
	 int headerLen = 0;
	 memset((void*)field, '\0', 1024);
	 memset((void*)number, '\0', 1024);
	 memset((void*)datatype, '\0', 1024);

	 parseHeader4(&headerLen, field, number, datatype);
	 if( myrank == 0) {
		 printf("field is %s, number is %s, 3rd field is %s, headerLen is %d\n", field, number, datatype, headerLen);
		 // skip header part
		 mfileCur += headerLen;
		 parsedTotalNumCells = (int) atoi(number);
		 // skip the cell data part as well, which is 5 int per cell
		 // TODO: check if it's valid integer (i.e. atoi does not handle long long
		 mfileCur += parsedTotalNumCells * sizeof(int) *5;
		 mfileCur ++; // skip extra '\n'
	 }
	 // rank 0 parsed the header of cell info, now bcast to everybody else
	 MPI_Bcast(&mfileCur, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);

	 // following code reads and skip cell type part
	 memset((void*)field, '\0', 1024);
	 memset((void*)number, '\0', 1024);
	 memset((void*)datatype, '\0', 1024);
	 parseHeader4short(&headerLen, field, number);
	 if( myrank == 0) {
		 printf("CELL type: field is %s, number is %s, 3rd field is %s, headerLen is %d\n",
				 field, number, datatype, headerLen);
		 // skip header part
		 mfileCur += headerLen;
		 parsedTotalNumCells = (int) atoi(number);
		 // skip the cell type data part as well, which is 1 int per cell
		 // TODO: check if it's valid integer (i.e. atoi does not handle long long
		 mfileCur += parsedTotalNumCells * sizeof(int) *1;
		 mfileCur ++; // skip extra '\n'
	 }
	 // rank 0 parsed the header of cell info, now bcast to everybody else
	 MPI_Bcast(&mfileCur, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);

	 // following code reads and skip POINT_DATA header
	 memset((void*)field, '\0', 1024);
	 memset((void*)number, '\0', 1024);
	 memset((void*)datatype, '\0', 1024);
	 parseHeader4short(&headerLen, field, number);
	 if( myrank == 0) {
		 printf("POINT_DATA: field is %s, number is %s, 3rd field is %s, headerLen is %d\n",
				 field, number, datatype, headerLen);
		 // skip header part
		 mfileCur += headerLen;
		 // it does not skip extra '\n' like previously, because POINT_DATA is a
		 // pure string and it does not have a data part followed by newline, like
		 // previous cases
		 //mfileCur ++; // skip extra '\n'
	 }
	 // rank 0 parsed the header of cell info, now bcast to everybody else
	 MPI_Bcast(&mfileCur, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
}

/**
 * This function reads 3D cell and cell type info
 *
 * @param:	numElems	#elements
 * @param:	numCells	#cells
 * @param:	numNodes	#nodes
 *
 */
#ifdef UPCASE
void READ3DCELLS4( int *eConnect, int *numElems, int *numCells, int *numNodes)
#elif  IBM
void read3dcells4( int *eConnect, int *numElems, int *numCells, int *numNodes)
#else
void read3dcells4_( int *eConnect, int *numElems, int *numCells, int *numNodes)
#endif
{
   int conn[9];
   int conn_new[9];
   int i, j;
   int elemType=12;

//	 long long totalNumCells = 0;
	 long long parsedTotalNumCells = 0;

	 char field[1024], number[1024], datatype[1024];
	 int headerLen = 0;
	 memset((void*)field, '\0', 1024);
	 memset((void*)number, '\0', 1024);
	 memset((void*)datatype, '\0', 1024);

	 parseHeader4(&headerLen, field, number, datatype);
	 if( myrank == 0) {
		 printf("field is %s, number is %s, 3rd field is %s, headerLen is %d\n", field, number, datatype, headerLen);
		 // skip header part
		 mfileCur += headerLen;
		 parsedTotalNumCells = (int) atoi(number);
		 // skip the cell data part as well, which is 9 int per cell
		 // TODO: check if it's valid integer (i.e. atoi does not handle long long
		 mfileCur += parsedTotalNumCells * sizeof(int) *9;
		 mfileCur ++; // skip extra '\n'
	 }
	 // rank 0 parsed the header of cell info, now bcast to everybody else
	 MPI_Bcast(&mfileCur, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);

	 // following code reads and skip cell type part
	 memset((void*)field, '\0', 1024);
	 memset((void*)number, '\0', 1024);
	 memset((void*)datatype, '\0', 1024);
	 parseHeader4short(&headerLen, field, number);
	 if( myrank == 0) {
		 printf("CELL type: field is %s, number is %s, 3rd field is %s, headerLen is %d\n",
				 field, number, datatype, headerLen);
		 // skip header part
		 mfileCur += headerLen;
		 parsedTotalNumCells = (int) atoi(number);
		 // skip the cell type data part as well, which is 1 int per cell
		 // TODO: check if it's valid integer (i.e. atoi does not handle long long
		 mfileCur += parsedTotalNumCells * sizeof(int) *1;
		 mfileCur ++; // skip extra '\n'
	 }
	 // rank 0 parsed the header of cell info, now bcast to everybody else
	 MPI_Bcast(&mfileCur, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);

	 // following code reads and skip POINT_DATA header
	 memset((void*)field, '\0', 1024);
	 memset((void*)number, '\0', 1024);
	 memset((void*)datatype, '\0', 1024);
	 parseHeader4short(&headerLen, field, number);
	 if( myrank == 0) {
		 printf("POINT_DATA: field is %s, number is %s, 3rd field is %s, headerLen is %d\n",
				 field, number, datatype, headerLen);
		 // skip header part
		 mfileCur += headerLen;
		 // it does not skip extra '\n' like previously, because POINT_DATA is a
		 // pure string and it does not have a data part followed by newline, like
		 // previous cases
		 //mfileCur ++; // skip extra '\n'
	 }
	 // rank 0 parsed the header of cell info, now bcast to everybody else
	 MPI_Bcast(&mfileCur, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
}


#ifdef UPCASE
void READFIELD4(int *fldid, double *vals, int *numNodes)
#elif  IBM
void readfield4(int *fldid, double *vals, int *numNodes)
#else
void readfield4_(int *fldid, double *vals, int *numNodes)
#endif
{
   int   i, j  ;
   char  fldname[100];

	 char field[1024], number[1024], datatype[1024];
	 int headerLen = 0;
	 memset((void*)field, '\0', 1024);
	 memset((void*)number, '\0', 1024);
	 memset((void*)datatype, '\0', 1024);

	 parseHeader4(&headerLen, field, number, datatype);
	 if( myrank == 0) {
		 printf("field1 is %s, field2 is %s, 3rd field is %s, headerLen is %d\n",
				 field, number, datatype, headerLen);
		 // skip header part
		 mfileCur += headerLen;
	 }
	 // rank 0 parsed the header of cell info, now bcast to everybody else
	 MPI_Bcast(&mfileCur, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);

	 long long llNumNodes = *numNodes;
	 long long totalNumNodes = 0;
	 long long myOffset;
	 MPI_Scan(&llNumNodes, &myOffset, 1, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);
	 myOffset -= (*numNodes);
	 myOffset *= sizeof(float) *3;

	 // get total numNodes
	 MPI_Allreduce( &llNumNodes, &totalNumNodes, 1, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);

	 float* fldval = (float*) malloc (sizeof(float) * 3* (*numNodes));
	 int ierr = MPI_File_read_at_all( mfile, mfileCur + myOffset, (void*)fldval,
			 (*numNodes)* sizeof(float) *3, MPI_CHAR, MPI_STATUSES_IGNORE);
	 if( ierr != 0) {
		 printf("error: MPI_File_read_at_all failed in readfield4()\n");
		 exit(4);
	 }
	 //update file pointer position
	 mfileCur += totalNumNodes * sizeof(float) * 3;
	 // +2 because there is an extra ' \n' at end of data block
	 mfileCur += 2;

	 for( i = 0; i < *numNodes; i++) {
		 swap_float_byte( &fldval[3*i+0]);
		 swap_float_byte( &fldval[3*i+1]);
		 swap_float_byte( &fldval[3*i+2]);
		 //printf("read field data: %E %E %E\n", fldval[3*i+0], fldval[3*i+1], fldval[3*i+2]);

		 vals[3*i+0] = fldval[3*i+0];
		 vals[3*i+1] = fldval[3*i+1];
		 vals[3*i+2] = fldval[3*i+2];
	 }
}
