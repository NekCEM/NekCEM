#ifndef H_IO_UTIL
#define H_IO_UTIL
#include <dirent.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/param.h>  // MAXPATHLEN definition

#define ONE_MILLION 1048576
#define kMaxPathLen MAXPATHLEN

extern int DEBUG_FLAG;
extern int IOTIMER_FLAG;
extern int IOTRACE_FLAG;

extern char filename[kMaxPathLen];

void getfieldname_( int i, char *name );
void getfilename_old(int *id, int *nid);
void adjust_endian();
int swap_int_byte(int *n);
int swap_float_byte(float *n);
int swap_long_long_byte(long long *n);

#endif

