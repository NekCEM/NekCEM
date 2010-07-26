#ifndef H_VTKBIN
#define H_VTKBIN
#include "rdtsc.h"

#define ONE_MILLION 1048576
#define BGP_FREQ (850*1000000)
#define V8_FREQ (3000*1000000)

extern char filename[100];
extern char mFilename[100];
extern char rbFilename[100];
extern char rbasciiFilename[100];

extern char* mfBuffer;
extern long long mfBufferCur ;
extern long long mfileCur ;
extern long long fieldSizeSum ;

extern int myrank, mysize;
extern int DEBUG_FLAG;

extern long long start_time;
extern long long end_time;
extern double overall_time;

void getfieldname_( int i, char *name );
void getfilename_(int *id, int *nid );
void adjust_endian();
int swap_int_byte(int *n);
int swap_float_byte(float *n);
int swap_long_long_byte(long long *n);


#endif

