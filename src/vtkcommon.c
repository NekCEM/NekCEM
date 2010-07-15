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

/**************************/

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

void getfilename_(int *id, int *nid )
{
   char ext0[100];
   char ext1[100];

/*printf( "\n  nid:: %d\n", *nid);*/
	memset((void*)filename, 0, 100);
	memset((void*)mFilename, 0, 100);
	memset((void*)rbFilename, 0, 100);

       strcpy( filename, "./vtk/binary-NN-p");
       sprintf( ext0, "%.6d-t", *nid);
       strcat( filename, ext0);
       sprintf( ext1, "%.5d", *id);
       strcat( filename, ext1);
       strcat( filename, ".vtk");

	strcpy( mFilename, "./vtk/mpi-binary-N1-t");
	sprintf( ext1, "%.5d", *id);
	strcat( mFilename, ext1);
	strcat( mFilename, ".vtk");

	strcpy (rbFilename, "./vtk/mpi-binary-NM1-t");
	sprintf( ext1, "%.5d", *id);
	strcat( rbFilename, ext1);
	strcat( rbFilename, ".vtk");

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
  /*
  printf("I'm little endian\n");
  */
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

