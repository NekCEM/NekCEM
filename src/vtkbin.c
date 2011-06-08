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
FILE *fp = NULL;

MPI_File mfile;

char filename[100];
char mFilename[100];

char* mfBuffer ;
long long mfileCur = 0, mfBufferCur = 0;
long long fieldSizeSum = 0;

//int myrank;



#ifdef UPCASE
void OPENFILE(  int *id, int *nid)
#elif  IBM
void openfile(  int *id, int *nid)
#else
void openfile_(  int *id, int *nid)
#endif
{
   getfilename_(id,nid, 3);
   fp = fopen(filename,  "w");
	 assert(fp);
}


#ifdef UPCASE
void CLOSEFILE()
#elif  IBM
void closefile()
#else
void closefile_()
#endif
{
   fclose(fp);
}

#ifdef UPCASE
void WRITEHEADER()
#elif  IBM
void writeheader()
#else
void writeheader_()
#endif
{
   int i ;/* np = 10;*/
   float xyz[3];
   assert( fp );

   /*printf("# vtk DataFile Version 2.0 \n"); */
   fprintf(fp, "# vtk DataFile Version 2.0 \n");
   fprintf(fp, "Electromagnetic Field  \n");
   fprintf(fp, "BINARY \n");
   fprintf(fp, "DATASET UNSTRUCTURED_GRID \n");
}

#ifdef UPCASE
void WRITENODES(double *xyzCoords, int *numNodes)
#elif  IBM
void writenodes(double *xyzCoords, int *numNodes)
#else
void writenodes_(double *xyzCoords, int *numNodes)
#endif
{
   float coord[3];
   int   i, j;
   fprintf(fp, "POINTS  %d ", *numNodes );
   fprintf(fp, " float  \n");
   for( i = 0; i < *numNodes; i++) {
       coord[0] = (float)xyzCoords[3*i+0];
       coord[1] = (float)xyzCoords[3*i+1];
       coord[2] = (float)xyzCoords[3*i+2];
       swap_float_byte( &coord[0] );
       swap_float_byte( &coord[1] );
       swap_float_byte( &coord[2] );
       fwrite(coord, sizeof(float), 3, fp);
   }
   fprintf( fp, " \n");
}

#ifdef UPCASE
void WRITE2DCELLS( int *eConnect, int *numElems, int *numCells, int *numNodes)
#elif  IBM
void write2dcells( int *eConnect, int *numElems, int *numCells, int *numNodes)
#else
void write2dcells_( int *eConnect, int *numElems, int *numCells, int *numNodes)
#endif
{
   int conn[5];
   int i, j;
   int elemType=9;

   fprintf( fp, "CELLS %d  %d \n", *numCells, 5*(*numCells));

   for (i = 0; i < *numCells; i++) {
        conn[0] = 4;
        conn[1] = eConnect[4*i+0];
        conn[2] = eConnect[4*i+1];
        conn[3] = eConnect[4*i+2];
        conn[4] = eConnect[4*i+3];
        for( j = 0; j < 5; j++) swap_int_byte( &conn[j] );
        fwrite(conn, sizeof(int), 5, fp);
   }
   fprintf( fp, "\n");
   fprintf( fp, "CELL_TYPES %d \n", *numCells);

   swap_int_byte(&elemType);

   for( i = 0; i < *numCells; i++)
    fwrite(&elemType,  sizeof(int), 1, fp);

   fprintf( fp, "\n");
   fprintf( fp, "POINT_DATA %d \n", *numNodes);
}

#ifdef UPCASE
void WRITE3DCELLS( int *eConnect, int *numElems, int *numCells, int *numNodes)
#elif  IBM
void write3dcells( int *eConnect, int *numElems, int *numCells, int *numNodes)
#else
void write3dcells_( int *eConnect, int *numElems, int *numCells, int *numNodes)
#endif
{
   int conn[9];
   int i, j;
   int elemType=12;

   fprintf( fp, "CELLS %d  %d \n", *numCells, 9*(*numCells) );

   for (i = 0; i < *numCells; i++) {
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
   }
   fprintf( fp, "\n");
   fprintf( fp, "CELL_TYPES %d \n", *numCells );

   swap_int_byte(&elemType);

   for (i = 0; i < *numCells; i++)
     fwrite(&elemType,  sizeof(int), 1, fp);

   fprintf( fp, "\n");
   fprintf( fp, "POINT_DATA %d \n", *numNodes );
}

#ifdef UPCASE
void WRITEFIELD(int *fldid, double *vals, int *numNodes)
#elif  IBM
void writefield(int *fldid, double *vals, int *numNodes)
#else
void writefield_(int *fldid, double *vals, int *numNodes)
#endif
{
   float fldval[3];
   int   i, j  ;
   char  fldname[100];

   getfieldname_(*fldid, fldname);

   fprintf( fp, "VECTORS %s ", fldname);
   fprintf( fp, " float \n");

   for (i = 0; i < *numNodes; i++) {

        fldval[0] = (float)vals[3*i+0];
        fldval[1] = (float)vals[3*i+1];
        fldval[2] = (float)vals[3*i+2];
        swap_float_byte( &fldval[0]);
        swap_float_byte( &fldval[1]);
        swap_float_byte( &fldval[2]);
        fwrite(fldval, sizeof(float), 3, fp);
   }
   fprintf(fp, " \n");
}


/*
 *  MPI-IO format (param(102)=4...) starts here
 */

#ifdef UPCASE
void OPENFILE4(  int *id, int *nid)
#elif  IBM
void openfile4(  int *id, int *nid)
#else
void openfile4_(  int *id, int *nid)
#endif
{
	getfilename_(id,nid, 4);

	/* parallel here*/

	int rc;
	rc = MPI_File_open(MPI_COMM_WORLD, mFilename, MPI_MODE_CREATE | MPI_MODE_RDWR , MPI_INFO_NULL, &mfile);
	if(rc){
		printf("Unable to create shared file %s in openfile\n", mFilename);
		fflush(stdout);
	}

	mfBuffer = (char*) malloc( sizeof( char) * 4 * ONE_MILLION);
	mfBufferCur = 0;

}

#ifdef UPCASE
void CLOSEFILE4()
#elif  IBM
void closefile4()
#else
void closefile4_()
#endif
{
   MPI_File_close( & mfile );
}


#ifdef UPCASE
void WRITEHEADER4(int* istep, int* idumpno, double* time, double* dt)
#elif  IBM
void writeheader4(int* istep, int* idumpno, double* time, double* dt)
#else
void writeheader4_(int* istep, int* idumpno, double* time, double* dt)
#endif
{
   int i ;/* np = 10;*/
   float xyz[3];

	 /*put header into sHeader string */
   char* sHeader = (char*)malloc(1024 * sizeof(char));
   memset((void*)sHeader, '\0', 1024);
   sprintf(sHeader, "# vtk DataFile Version 3.0 \n");
//   sprintf(sHeader+strlen(sHeader), "Electromagnetic Field  \n");
   sprintf(sHeader+strlen(sHeader), "%d %d %E %E Electromagnetic Field  \n",
					 *istep, *idumpno, *time, *dt);
   sprintf(sHeader+strlen(sHeader),  "BINARY \n");
   sprintf(sHeader+strlen(sHeader), "DATASET UNSTRUCTURED_GRID \n");


   mfileCur = 0; //reset file position for new file
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank );
   if(myrank == 0) //HARDCODED for now, should be first local_rank
   {
   memcpy(mfBuffer, sHeader, strlen(sHeader));
   mfBufferCur = strlen(sHeader);
   }

   free(sHeader);
//printf("writeheader4() done\n");

}


#ifdef UPCASE
void WRITENODES4(double *xyzCoords, int *numNodes)
#elif  IBM
void writenodes4(double *xyzCoords, int *numNodes)
#else
void writenodes4_(double *xyzCoords, int *numNodes)
#endif
{
   float coord[3];
   int   i, j;
/*
   fprintf(fp, "POINTS  %d ", *numNodes );
   fprintf(fp, " float  \n");
*/
	//has to change *numNodes here
	//numNodes would be aggregated
	//xyzCoords itself doesn't change (absolute value)
	int totalNumNodes = 0;
	MPI_Allreduce(numNodes, &totalNumNodes, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD);
	if( myrank == 0)
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
       swap_float_byte( &coord[0] );
       swap_float_byte( &coord[1] );
       swap_float_byte( &coord[2] );
//       fwrite(coord, sizeof(float), 3, fp);

	memcpy(&mfBuffer[mfBufferCur], coord, sizeof(float)*3);
	mfBufferCur += sizeof(float) *3;
   }
//   fprintf( fp, " \n");


//	mfBuffer[mfBufferCur++] = '\n';

	//my_data_offset is local offset for big chunk of field/cell data
	//mfileCur is the glocal file offset shared with every proc
	long long my_data_offset = 0;
	MPI_Status write_status;
	MPI_Scan(&mfBufferCur, &my_data_offset, 1, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&mfBufferCur, &fieldSizeSum,  1, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);

	my_data_offset += mfileCur;
	my_data_offset -= mfBufferCur;
	MPI_File_write_at_all_begin(mfile, my_data_offset, mfBuffer, mfBufferCur, MPI_CHAR);
	MPI_File_write_at_all_end(mfile, mfBuffer, &write_status);

	mfileCur += fieldSizeSum;
//	MPI_File_close( & mfile );

	//clean up mfBuffer/mfBufferCur
	mfBufferCur = 0;
	//simply adding "\n", following original format
	if( myrank == 0)
	{
	mfBuffer[mfBufferCur] = '\n';
	mfBufferCur++;
	}
//printf("writenodes4() done\n");

}

#ifdef UPCASE
void WRITE2DCELLS4( int *eConnect, int *numElems, int *numCells, int *numNodes)
#elif  IBM
void write2dcells4( int *eConnect, int *numElems, int *numCells, int *numNodes)
#else
void write2dcells4_( int *eConnect, int *numElems, int *numCells, int *numNodes)
#endif
{
   int conn[5];
   int conn_new[5];
   int i, j;
   int elemType=9;

//   fprintf( fp, "CELLS %d  %d \n", *numCells, 5*(*numCells));

//cell number would be aggregated here
//following conn number would add an offset - myeConnOffset
        int totalNumCells = 0;
        MPI_Allreduce( numCells, &totalNumCells, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD);

        int myeConnOffset = 0;
        MPI_Scan( numNodes, &myeConnOffset, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD);
        myeConnOffset -= *numNodes;

        int totalNumNodes = 0;
        MPI_Allreduce(numNodes, &totalNumNodes, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD);

        if( myrank == 0)
        {
        char* sHeader = (char*) malloc (1024 * sizeof(char));
        memset((void*)sHeader, '\0', 1024);
        sprintf(sHeader, "CELLS %d  %d \n", totalNumCells, 5*(totalNumCells) );

        memcpy(&mfBuffer[mfBufferCur], sHeader, strlen(sHeader));
        mfBufferCur += strlen(sHeader);
        free(sHeader);
        }


   for (i = 0; i < *numCells; i++) {
/*
        conn[0] = 4;
        conn[1] = eConnect[4*i+0];
        conn[2] = eConnect[4*i+1];
        conn[3] = eConnect[4*i+2];
        conn[4] = eConnect[4*i+3];
	for( j = 0; j < 5; j++) swap_int_byte( &conn[j] );
        fwrite(conn, sizeof(int), 5, fp);
*/
//mpi-io part
        conn_new[0] = 4;
        conn_new[1] = eConnect[4*i+0] + myeConnOffset;
        conn_new[2] = eConnect[4*i+1] + myeConnOffset;
        conn_new[3] = eConnect[4*i+2] + myeConnOffset;
        conn_new[4] = eConnect[4*i+3] + myeConnOffset;
        for( j = 0; j < 5; j++) swap_int_byte( &conn_new[j] );

        memcpy(&mfBuffer[mfBufferCur], conn_new, sizeof(int)*5);
        mfBufferCur += sizeof(int) * 5;

   }
//flush to disk

        long long my_data_offset = 0;
        MPI_Status write_status;
        MPI_Scan(&mfBufferCur, &my_data_offset, 1, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&mfBufferCur, &fieldSizeSum,  1, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);

        my_data_offset += mfileCur;
        my_data_offset -= mfBufferCur;
        MPI_File_write_at_all_begin(mfile, my_data_offset, mfBuffer, mfBufferCur, MPI_CHAR);
        MPI_File_write_at_all_end(mfile, mfBuffer, &write_status);

        mfileCur += fieldSizeSum;
        mfBufferCur = 0;
/*
   fprintf( fp, "\n");
   fprintf( fp, "CELL_TYPES %d \n", *numCells);
*/
//mpi-io part

	if( myrank == 0)
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
//    fwrite(&elemType,  sizeof(int), 1, fp);

//mpi-io
    memcpy(&mfBuffer[mfBufferCur], &elemType, sizeof(int));
    mfBufferCur += sizeof(int);
   }

	my_data_offset = 0;
        MPI_Scan(&mfBufferCur, &my_data_offset, 1, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&mfBufferCur, &fieldSizeSum,  1, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);

        my_data_offset += mfileCur;
        my_data_offset -= mfBufferCur;
        MPI_File_write_at_all_begin(mfile, my_data_offset, mfBuffer, mfBufferCur, MPI_CHAR);
        MPI_File_write_at_all_end(mfile, mfBuffer, &write_status);

        mfileCur += fieldSizeSum;
        mfBufferCur = 0;
/*
   fprintf( fp, "\n");
   fprintf( fp, "POINT_DATA %d \n", *numNodes);
*/

//mpi-io

	if( myrank == 0)
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
void WRITE3DCELLS4( int *eConnect, int *numElems, int *numCells, int *numNodes)
#elif  IBM
void write3dcells4( int *eConnect, int *numElems, int *numCells, int *numNodes)
#else
void write3dcells4_( int *eConnect, int *numElems, int *numCells, int *numNodes)
#endif
{
   int conn[9];
   int conn_new[9];
   int i, j;
   int elemType=12;
//   fprintf( fp, "CELLS %d  %d \n", *numCells, 9*(*numCells) );

//cell number would be aggregated here
//following conn number would add an offset - myeConnOffset
//
	int totalNumCells = 0;
	MPI_Allreduce( numCells, &totalNumCells, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD);

	int myeConnOffset = 0;
	MPI_Scan( numNodes, &myeConnOffset, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD);
	myeConnOffset -= *numNodes;

	int totalNumNodes = 0;
        MPI_Allreduce(numNodes, &totalNumNodes, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD);

	if( myrank == 0)
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
//mpi-io part
	conn_new[0] = 8;
        conn_new[1] = eConnect[8*i+0] + myeConnOffset;
        conn_new[2] = eConnect[8*i+1] + myeConnOffset;
        conn_new[3] = eConnect[8*i+2] + myeConnOffset;
        conn_new[4] = eConnect[8*i+3] + myeConnOffset;
        conn_new[5] = eConnect[8*i+4] + myeConnOffset;
        conn_new[6] = eConnect[8*i+5] + myeConnOffset;
        conn_new[7] = eConnect[8*i+6] + myeConnOffset;
        conn_new[8] = eConnect[8*i+7] + myeConnOffset;
	for( j = 0; j < 9; j++) swap_int_byte( &conn_new[j] );

	memcpy(&mfBuffer[mfBufferCur], conn_new, sizeof(int)*9);
	mfBufferCur += sizeof(int) * 9;
   }

	long long my_data_offset = 0;
        MPI_Status write_status;
        MPI_Scan(&mfBufferCur, &my_data_offset, 1, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&mfBufferCur, &fieldSizeSum,  1, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);

        my_data_offset += mfileCur;
        my_data_offset -= mfBufferCur;
        MPI_File_write_at_all_begin(mfile, my_data_offset, mfBuffer, mfBufferCur, MPI_CHAR);
        MPI_File_write_at_all_end(mfile, mfBuffer, &write_status);

        mfileCur += fieldSizeSum;
        mfBufferCur = 0;
/*
   fprintf( fp, "\n");
   fprintf( fp, "CELL_TYPES %d \n", *numCells );
*/
	if( myrank == 0)
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
//     fwrite(&elemType,  sizeof(int), 1, fp);

	memcpy(&mfBuffer[mfBufferCur], &elemType, sizeof(int));
	mfBufferCur += sizeof(int);
   }

	my_data_offset = 0;
        MPI_Scan(&mfBufferCur, &my_data_offset, 1, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&mfBufferCur, &fieldSizeSum,  1, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);

        my_data_offset += mfileCur;
        my_data_offset -= mfBufferCur;
        MPI_File_write_at_all_begin(mfile, my_data_offset, mfBuffer, mfBufferCur, MPI_CHAR);
        MPI_File_write_at_all_end(mfile, mfBuffer, &write_status);

        mfileCur += fieldSizeSum;
	mfBufferCur = 0;

/*
   fprintf( fp, "\n");
   fprintf( fp, "POINT_DATA %d \n", *numNodes );
*/
	if( myrank == 0)
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
void WRITEFIELD4(int *fldid, double *vals, int *numNodes)
#elif  IBM
void writefield4(int *fldid, double *vals, int *numNodes)
#else
void writefield4_(int *fldid, double *vals, int *numNodes)
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

	if( myrank == 0)
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
        swap_float_byte( &fldval[0]);
        swap_float_byte( &fldval[1]);
        swap_float_byte( &fldval[2]);
//        fwrite(fldval, sizeof(float), 3, fp);

	memcpy( &mfBuffer[mfBufferCur], fldval, sizeof(float)*3);
	mfBufferCur += sizeof(float) *3;
   }
	long long my_data_offset = 0;
        MPI_Status write_status;
        MPI_Scan(&mfBufferCur, &my_data_offset, 1, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&mfBufferCur, &fieldSizeSum,  1, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);

        my_data_offset += mfileCur;
        my_data_offset -= mfBufferCur;
        MPI_File_write_at_all_begin(mfile, my_data_offset, mfBuffer, mfBufferCur, MPI_CHAR);
        MPI_File_write_at_all_end(mfile, mfBuffer, &write_status);

        mfileCur += fieldSizeSum;
        mfBufferCur = 0;

//   fprintf(fp, " \n");

	//add this return symbol into mpi file...
	if( myrank == 0)
        {
        char* sHeader = (char*) malloc (1024 * sizeof(char));
        memset((void*)sHeader, '\0', 1024);
        sprintf(sHeader, " \n");

        memcpy(&mfBuffer[mfBufferCur], sHeader, strlen(sHeader));
        mfBufferCur += strlen(sHeader);
        free(sHeader);
        }
}

/*********************************************************/
// following is for element mapping file I/O
/*********************************************************/
#define MAPPING_FILE_NAME "./vtk/element_numbering.dat"
#define MAPPING_FILE_HEADER_SIZE 1024

/**
 * This function takes the global array and size, and write to a file
 * (wrapper write func)
 *
 * eltNum: glboal element numbering array
 * nelt:	 size of eleNum
 */
void file_write_element_numbering(  int *eltNum, int *nelt)
{
	// open binary mapping file
	FILE* fp = NULL;
	fp = fopen(MAPPING_FILE_NAME, "wb");
	if (fp == NULL) { printf("error: can't open file %s \n", MAPPING_FILE_NAME); exit(1);}

	char mapping_header_str[MAPPING_FILE_HEADER_SIZE];
	memset(mapping_header_str, '\0', MAPPING_FILE_HEADER_SIZE);

	int g_size;
	MPI_Comm_size(MPI_COMM_WORLD, &g_size);
	sprintf(mapping_header_str, "processor number: %d", g_size);

	int ret;
	ret = fwrite(mapping_header_str, sizeof(char), MAPPING_FILE_HEADER_SIZE, fp);
	assert(ret == MAPPING_FILE_HEADER_SIZE);
	ret = fwrite(eltNum, sizeof(int), (size_t) (*nelt), fp);
	assert(ret == (*nelt));

	fclose(fp);
}

/**
 * This function takes local element numbering array and size, assemble them,
 * and write to a single file in order of processor ranks
 *
 * local_elm: local element numbering array
 * nelt:			size of local_elm
 */
#ifdef UPCASE
void WRITE_ELEMENT_NUMBERING(  int *local_elm, int *nelt)
#elif  IBM
void write_element_numbering(  int *local_elm, int *nelt)
#else
void write_element_numbering_(  int *local_elm, int *nelt)
#endif
{
	int gsize;
	int root = 0;
	int *displs, i, *rcounts, *rbuf;

	int mysize, myrank;
	MPI_Comm_size(MPI_COMM_WORLD, &gsize);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank );
	displs = (int*) malloc(sizeof(int) * gsize);
	rcounts = (int*) malloc(sizeof(int) * gsize);
	assert(displs != NULL);
	assert(rcounts != NULL);
	// gather all sizes of local_elm array
	MPI_Gather(nelt, 1, MPI_INT, rcounts, 1, MPI_INT, root, MPI_COMM_WORLD);

	//transform local size array to displs
	displs[0] = 0;
	for (i = 1; i < gsize; i++) {
		displs[i] = displs[i-1] + rcounts[i-1];
	}

	// alloc recv buffer
	int total_nelt;
	MPI_Allreduce(nelt, &total_nelt, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	rbuf = (int*) malloc(sizeof(int) * total_nelt);
	assert(rbuf != NULL);

	MPI_Gatherv(local_elm, *nelt, MPI_INT, rbuf, rcounts, displs, MPI_INT, root, MPI_COMM_WORLD);

	if (myrank == root) {
		file_write_element_numbering(rbuf, &total_nelt);
	}

	MPI_Barrier(MPI_COMM_WORLD);
}

/**
 * This function takes size of array, and read them to the array from a file
 * (wrapper read func)
 *
 * eltNum: glboal element numbering array (output)
 * nelt:	 size of eleNum (input)
 */
void file_read_element_numbering(  int *eltNum, int *nelt)
{
	// open binary mapping file
	FILE* fp = NULL;
	fp = fopen(MAPPING_FILE_NAME, "rb");
	if (fp == NULL) { printf("error: can't open file %s \n", MAPPING_FILE_NAME); exit(1);}

	int ret;
	ret = fseek (fp, MAPPING_FILE_HEADER_SIZE, SEEK_SET);
	assert(ret == 0);

	ret = fread(eltNum, sizeof(int), (size_t) (*nelt), fp);
	assert(ret == (*nelt));

	fclose(fp);
}

/**
 * This function takes size of global element array, reads the element
 * numbering from a data file, and put it in global element numbering array
 *
 * eltNum:	global element numbering array (output)
 * nelt:		size of eltNum (input)
 */
#ifdef UPCASE
void READ_ELEMENT_NUMBERING(  int *eltNum, int *nelt)
#elif  IBM
void read_element_numbering(  int *eltNum, int *nelt)
#else
void read_element_numbering_(  int *eltNum, int *nelt)
#endif
{
	int root = 0;

	int mysize, myrank;
	MPI_Comm_size(MPI_COMM_WORLD, &mysize);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank );

	if( myrank == root ) {
		file_read_element_numbering(eltNum, nelt);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	// bcast to every other processor
	MPI_Bcast(eltNum, *nelt, MPI_INT, root, MPI_COMM_WORLD);
}




