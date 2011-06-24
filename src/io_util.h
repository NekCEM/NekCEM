#ifndef H_IO_UTIL
#define H_IO_UTIL
#include <dirent.h>
#include <unistd.h>
#include <sys/stat.h>

#define ONE_MILLION 1048576

extern int DEBUG_FLAG;
extern int IOTIMER_FLAG;
extern int IOTRACE_FLAG;

extern char filename[100];

void getfieldname_( int i, char *name );
void getfilename_old(int *id, int *nid);
void adjust_endian();
int swap_int_byte(int *n);
int swap_float_byte(float *n);
int swap_long_long_byte(long long *n);

#endif

