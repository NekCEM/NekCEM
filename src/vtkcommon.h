#ifndef H_VTKBIN
#define H_VTKBIN
#include <dirent.h>
#include <unistd.h>
#include <sys/stat.h>

#include "rdtsc.h"
#include <mpi.h>

#define ONE_MILLION 1048576
#define BGP_FREQ (850*1000000)
#define V8_FREQ (3000*1000000)
#define WRITERBUFFERSIZE (50 * ONE_MILLION)

#define kStrLocal "./vtk"
#define kStrFs0Misun "/intrepid-fs0/users/mmin/scratch/NEKCEM_vtk"
#define kStrFs0Fuji "/intrepid-fs1/users/fuji/scratch/NEKCEM_vtk"

//#define kOutputPath kStrLocal
#define kOutputPath kStrLocal
//#define kOutputPath kStrFs0Misun

extern char filename[100];
extern char mFilename[100];
extern char rbFilename[100];
extern char rbasciiFilename[100];
extern char rbnmmFilename[128];
extern char nmFilename[128];

extern char* mfBuffer;
extern long long mfBufferCur ;
extern long long mfileCur ;
extern long long fieldSizeSum ;

extern int myrank, mysize;
extern int mySpecies, localsize, groupRank, numGroups;
extern MPI_Comm localcomm;

extern int DEBUG_FLAG;
extern int IOTIMER_FLAG;
extern int IOTRACE_FLAG;

extern long long start_time;
extern long long end_time;
extern double overall_time;

void getfieldname_( int i, char *name );
void getfilename_(int *id, int *nid, int io_option);
void adjust_endian();
int swap_int_byte(int *n);
int swap_float_byte(float *n);
int swap_long_long_byte(long long *n);

#endif

