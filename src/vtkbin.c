/**
 * This file contains functions for param(81) =3, ie. the legacy C-POSIX I/O
 * library that generates 1PFPP. It inlcudes io_util or mpiio_util.
 *
 * These function can be called in either MPI or non-MPI environments.
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#ifdef MPI
#include "mpiio_util.h"
#else
#include "io_util.h"
#endif

FILE *fp1 = NULL;

char filename[kMaxPathLen];

#ifdef UPCASE
void OPENFILE(  int *id, int *nid)
#elif  IBM
void openfile(  int *id, int *nid)
#else
void openfile_(  int *id, int *nid)
#endif
{
	// add following ifdef because io_option 3 could run on either non-mpi or
	// mpi
#ifdef MPI
   getfilename_(id,nid, 3);
#else
	 getfilename_old(id, nid);
#endif
   fp1 = fopen(filename,  "w");
	 assert(fp1);
}

#ifdef UPCASE
void CLOSEFILE()
#elif  IBM
void closefile()
#else
void closefile_()
#endif
{
   fclose(fp1);
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
   assert( fp1 );

   /*printf("# vtk DataFile Version 2.0 \n"); */
   fprintf(fp1, "# vtk DataFile Version 2.0 \n");
   fprintf(fp1, "Electromagnetic Field  \n");
   fprintf(fp1, "BINARY \n");
   fprintf(fp1, "DATASET UNSTRUCTURED_GRID \n");
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
   fprintf(fp1, "POINTS  %d ", *numNodes );
   fprintf(fp1, " float  \n");
   for( i = 0; i < *numNodes; i++) {
       coord[0] = (float)xyzCoords[3*i+0];
       coord[1] = (float)xyzCoords[3*i+1];
       coord[2] = (float)xyzCoords[3*i+2];
       swap_float_byte( &coord[0] );
       swap_float_byte( &coord[1] );
       swap_float_byte( &coord[2] );
       fwrite(coord, sizeof(float), 3, fp1);
   }
   fprintf( fp1, " \n");
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

   fprintf( fp1, "CELLS %d  %d \n", *numCells, 5*(*numCells));

   for (i = 0; i < *numCells; i++) {
        conn[0] = 4;
        conn[1] = eConnect[4*i+0];
        conn[2] = eConnect[4*i+1];
        conn[3] = eConnect[4*i+2];
        conn[4] = eConnect[4*i+3];
        for( j = 0; j < 5; j++) swap_int_byte( &conn[j] );
        fwrite(conn, sizeof(int), 5, fp1);
   }
   fprintf( fp1, "\n");
   fprintf( fp1, "CELL_TYPES %d \n", *numCells);

   swap_int_byte(&elemType);

   for( i = 0; i < *numCells; i++)
    fwrite(&elemType,  sizeof(int), 1, fp1);

   fprintf( fp1, "\n");
   fprintf( fp1, "POINT_DATA %d \n", *numNodes);
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

   fprintf( fp1, "CELLS %d  %d \n", *numCells, 9*(*numCells) );

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
        fwrite(conn, sizeof(int), 9, fp1);
   }
   fprintf( fp1, "\n");
   fprintf( fp1, "CELL_TYPES %d \n", *numCells );

   swap_int_byte(&elemType);

   for (i = 0; i < *numCells; i++)
     fwrite(&elemType,  sizeof(int), 1, fp1);

   fprintf( fp1, "\n");
   fprintf( fp1, "POINT_DATA %d \n", *numNodes );
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

   fprintf( fp1, "VECTORS %s ", fldname);
   fprintf( fp1, " float \n");

   for (i = 0; i < *numNodes; i++) {

        fldval[0] = (float)vals[3*i+0];
        fldval[1] = (float)vals[3*i+1];
        fldval[2] = (float)vals[3*i+2];
        swap_float_byte( &fldval[0]);
        swap_float_byte( &fldval[1]);
        swap_float_byte( &fldval[2]);
        fwrite(fldval, sizeof(float), 3, fp1);
   }
   fprintf(fp1, " \n");
}


#ifndef MPI

#ifdef UPCASE
void OPENFILE4(  int *id, int *nid);
#elif  IBM
void openfile4(  int *id, int *nid);
#else
void openfile4_(  int *id, int *nid){}
#endif


#ifdef UPCASE
void CLOSEFILE4();
#elif  IBM
void closefile4();
#else
void closefile4_(){};
#endif

#ifdef UPCASE
void WRITEHEADER4(int* istep, int* idumpno, double* time, double* dt);
#elif  IBM
void writeheader4(int* istep, int* idumpno, double* time, double* dt);
#else
void writeheader4_(int* istep, int* idumpno, double* time, double* dt){};
#endif


#ifdef UPCASE
void WRITENODES4(double *xyzCoords, int *numNodes);
#elif  IBM
void writenodes4(double *xyzCoords, int *numNodes);
#else
void writenodes4_(double *xyzCoords, int *numNodes){};
#endif


#ifdef UPCASE
void WRITE2DCELLS4( int *eConnect, int *numElems, int *numCells, int *numNodes);
#elif  IBM
void write2dcells4( int *eConnect, int *numElems, int *numCells, int *numNodes);
#else
void write2dcells4_( int *eConnect, int *numElems, int *numCells, int *numNodes){};
#endif


#ifdef UPCASE
void WRITE2DCELLS4_SWAP( int *eConnect, int *numElems, int *numCells, int *numNodes);
#elif  IBM
void write2dcells4_swap( int *eConnect, int *numElems, int *numCells, int *numNodes);
#else
void write2dcells4_swap_( int *eConnect, int *numElems, int *numCells, int *numNodes){};
#endif

#ifdef UPCASE
void WRITE3DCELLS4( int *eConnect, int *numElems, int *numCells, int *numNodes);
#elif  IBM
void write3dcells4( int *eConnect, int *numElems, int *numCells, int *numNodes);
#else
void write3dcells4_( int *eConnect, int *numElems, int *numCells, int *numNodes){};
#endif

#ifdef UPCASE
void WRITECELLS4( int *eConnect, int *numElems, int *numCells, int *numNodes);
#elif  IBM
void writecells4( int *eConnect, int *numElems, int *numCells, int *numNodes);
#else
void writecells4_( int *eConnect, int *numElems, int *numCells, int *numNodes){};
#endif

#ifdef UPCASE
void WRITE3DCELLS4_SWAP( int *eConnect, int *numElems, int *numCells, int *numNodes);
#elif  IBM
void write3dcells4_swap( int *eConnect, int *numElems, int *numCells, int *numNodes);
#else
void write3dcells4_swap_( int *eConnect, int *numElems, int *numCells, int *numNodes){};
#endif

#ifdef UPCASE
void WRITEFIELD4(int *fldid, double *vals, int *numNodes);
#elif  IBM
void writefield4(int *fldid, double *vals, int *numNodes);
#else
void writefield4_(int *fldid, double *vals, int *numNodes){};
#endif

#ifdef UPCASE
void WRITEFIELD4_DOUBLE(int *fldid, double *vals, int *numNodes);
#elif  IBM
void writefield4_double(int *fldid, double *vals, int *numNodes);
#else
void writefield4_double_(int *fldid, double *vals, int *numNodes){};
#endif

#ifdef UPCASE
void INITRBIO(int *numgroups)
#elif  IBM
void initrbio(int *numgroups)
#else
void initrbio_(int *numgroups)
#endif
{}

#ifdef UPCASE
void OPENFILE6(  int *id, int *nid)
#elif  IBM
void openfile6(  int *id, int *nid)
#else
void openfile6_(  int *id, int *nid)
#endif
{}

#ifdef UPCASE
void CLOSEFILE6()
#elif  IBM
void closefile6()
#else
void closefile6_()
#endif
{}

#ifdef UPCASE
void OPENFILE_RESTART(  int *id, int *nid)
#elif  IBM
void openfile_restart(  int *id, int *nid)
#else
void openfile_restart_(  int *id, int *nid)
#endif
{}

#ifdef UPCASE
void CLOSEFILE_RESTART()
#elif  IBM
void closefile_restart()
#else
void closefile_restart_()
#endif
{}

#ifdef UPCASE
void WRITEHEADER6()
#elif  IBM
void writeheader6()
#else
void writeheader6_()
#endif
{}

#ifdef UPCASE
void WRITENODES6(double *xyzCoords, int *numNodes)
#elif  IBM
void writenodes6(double *xyzCoords, int *numNodes)
#else
void writenodes6_(double *xyzCoords, int *numNodes)
#endif
{}

#ifdef UPCASE
void WRITE2DCELLS6( int *eConnect, int *numElems, int *numCells, int *numNodes)
#elif  IBM
void write2dcells6( int *eConnect, int *numElems, int *numCells, int *numNodes)
#else
void write2dcells6_( int *eConnect, int *numElems, int *numCells, int *numNodes)
#endif
{}

#ifdef UPCASE
void WRITE3DCELLS6( int *eConnect, int *numElems, int *numCells, int *numNodes)
#elif  IBM
void write3dcells6( int *eConnect, int *numElems, int *numCells, int *numNodes)
#else
void write3dcells6_( int *eConnect, int *numElems, int *numCells, int *numNodes)
#endif
{}

#ifdef UPCASE
void WRITEFIELD6(int *fldid, double *vals, int *numNodes)
#elif  IBM
void writefield6(int *fldid, double *vals, int *numNodes)
#else
void writefield6_(int *fldid, double *vals, int *numNodes)
#endif
{}

#ifdef UPCASE
void SET_ASCII_TRUE(int *numgroups)
#elif  IBM
void set_ascii_true(int *numgroups)
#else
void set_ascii_true_(int *numgroups)
#endif
{}

#ifdef UPCASE
void SET_ASCII_NMM(int *numgroups)
#elif  IBM
void set_ascii_nmm(int *numgroups)
#else
void set_ascii_nmm_(int *numgroups)
#endif
{}

#ifdef UPCASE
void SET_ASCII_NM(int *numgroups)
#elif  IBM
void set_ascii_nm(int *numgroups)
#else
void set_ascii_nm_(int *numgroups)
#endif
{}

#ifdef UPCASE
void PVTK_NMM(  int *id)
#elif  IBM
void pvtk_nmm(  int *id)
#else
void pvtk_nmm_(  int *id)
#endif
{}

#ifdef UPCASE
void STARTTIMING()
#elif  IBM
void starttiming()
#else
void starttiming_()
#endif
{}


#ifdef UPCASE
void ENDTIMING()
#elif  IBM
void endtiming()
#else
void endtiming_()
#endif
{}

#ifdef UPCASE
void PASS_IO_PARAMS(int *param1, int* param2)
#elif  IBM
void pass_io_params(int *param1, int* param2)
#else
void pass_io_params_(int *param1, int* param2)
#endif
{}

#ifdef UPCASE
void PRINTIO(int *fp1aram, int* piostep)
#elif  IBM
void printio(int *fp1aram, int* piostep)
#else
void printio_(int *fp1aram, int* piostep)
#endif
{}

#ifdef UPCASE
void WRITEIOTRACE(int *fp1aram, int* piostep)
#elif  IBM
void writeiotrace(int *fp1aram, int* piostep)
#else
void writeiotrace_(int *fp1aram, int* piostep)
#endif
{}

#ifdef UPCASE
void WRITECOMPUTETRACE(int *fp1aram, int* pnf, int* pcompstep, double* pdtime, double* pcpu_t)
#elif  IBM
void writecomputetrace(int *fp1aram, int* pnf, int* pcompstep, double* pdtime, double* pcpu_t)
#else
void writecomputetrace_(int *fp1aram, int* pnf, int* pcompstep, double* pdtime, double* pcpu_t)
#endif
{}

#ifdef UPCASE
void FREE_RBIO_BUFFER ()
#elif  IBM
void free_rbio_buffer ()
#else
void free_rbio_buffer_()
#endif
{}

#endif
