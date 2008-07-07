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


FILE *fp = NULL;
char filename[100];

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

void getfilename_(int *id, int *nid )
{
   char ext0[100];
   char ext1[100];

       /*printf( "\n  nid:: %d\n", *nid);*/
       strcpy( filename, "./vtk/em-p");
       sprintf( ext0, "%.6d-t", *nid);
       strcat( filename, ext0);
       sprintf( ext1, "%.5d", *id);
       strcat( filename, ext1);
       strcat( filename, ".vtk");
}

int swap_int_byte(int *n)
{
  unsigned char *cptr,tmp;

  cptr = (unsigned char *)n;
  tmp = cptr[0];
  cptr[0] = cptr[3];
  cptr[3] = tmp;
  tmp = cptr[1];
  cptr[1] = cptr[2];
  cptr[2] = tmp;

  return 0;
}

int swap_float_byte(float *n)
{
  unsigned char *cptr,tmp;

  cptr = (unsigned char *)n;
  tmp  = cptr[0];
  cptr[0] = cptr[3];
  cptr[3] = tmp    ;
  tmp     = cptr[1];
  cptr[1] = cptr[2];
  cptr[2] = tmp    ;
  return 0;
}

#ifdef UPCASE
void OPENFILE(  int *id, int *nid)
#elif  IBM
void openfile(  int *id, int *nid)
#else
void openfile_(  int *id, int *nid)
#endif
{
   getfilename_(id,nid);
   fp = fopen(filename,  "w"); assert(fp);
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
                   


