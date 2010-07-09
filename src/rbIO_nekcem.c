#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include <mpi.h>
#include "vtkbin.h"



MPI_File mfile;
char mFilename[100];

char mfBuffer[ 4* ONE_MILLION];
long long mfileCur = 0, mfBufferCur = 0;
long long fieldSizeSum = 0;

int myrank;

/*
 *  MPI-IO format (param(102)=4...) starts here
 */

#ifdef UPCASE
void OPENFILE6(  int *id, int *nid)
#elif  IBM
void openfile6(  int *id, int *nid)
#else
void openfile6_(  int *id, int *nid)
#endif
{
   getfilename_(id,nid);

/* parallel here*/

        int rc;
        rc = MPI_File_open(MPI_COMM_WORLD, mFilename, MPI_MODE_CREATE | MPI_MODE_RDWR , MPI_INFO_NULL, &mfile);
        if(rc){
            printf("Unable to create shared file %s in openfile\n", mFilename);
            fflush(stdout);
          }
        mfBufferCur = 0;

}
