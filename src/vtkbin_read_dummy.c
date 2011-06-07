#include <stdio.h>
#include <stdlib.h>
#include <assert.h>


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
void OPENFILE_READ4(int *id,  int *nid)
#elif  IBM
void openfile_read4(int *id,  int *nid)
#else
void openfile_read4_(int *id,  int *nid){}
#endif

#ifdef UPCASE
void CLOSEFILE_READ4()
#elif  IBM
void closefile_read4()
#else
void closefile_read4_(){}
#endif


#ifdef UPCASE
void READHEADER4(int* irsttemp, int* idump)
#elif  IBM
void readheader4(int* irsttemp, int* idump)
#else
void readheader4_(int* irsttemp, int* idump){}
#endif

/**
 * This function parse the data header that has only two parts, e.g. cell
 * type, points, etc.
 *
 */
void parseHeader4short(int* len, char* field, char* number) {
}

/**
 * This function parse the data header that has three parts, e.g. nodes
 *
 */
void parseHeader4(int* len, char* field, char* number, char* datatype) {
}

#ifdef UPCASE
void READNODES4(double *xyzCoords, int *numNodes)
#elif  IBM
void readnodes4(double *xyzCoords, int *numNodes)
#else
void readnodes4_(double *xyzCoords, int *numNodes){}
#endif


#ifdef UPCASE
void READ2DCELLS4( int *eConnect, int *numElems, int *numCells, int *numNodes)
#elif  IBM
void read2dcells4( int *eConnect, int *numElems, int *numCells, int *numNodes)
#else
void read2dcells4_( int *eConnect, int *numElems, int *numCells, int *numNodes){}
#endif

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
void read3dcells4_( int *eConnect, int *numElems, int *numCells, int *numNodes){}
#endif

#ifdef UPCASE
void READFIELD4(int *fldid, double *vals, int *numNodes)
#elif  IBM
void readfield4(int *fldid, double *vals, int *numNodes)
#else
void readfield4_(int *fldid, double *vals, int *numNodes){}
#endif
