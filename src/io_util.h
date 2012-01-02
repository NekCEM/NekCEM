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

// This is for writer's buffer size in rbIO
//#define WRITERBUFFERSIZE (50 * ONE_MILLION)
extern long long WRITERBUFFERSIZE;

// Specify output path
#define kStrLocal "./vtk"
#define kStrFs0Misun "/intrepid-fs0/users/mmin/scratch/NEKCEM_vtk"
#define kStrFs0Fuji "/intrepid-fs0/users/fuji/scratch/NEKCEM_vtk"
#define kStrTitanJing "/tmp/work/jingfu/NEKCEM_vtk"

extern char* kOutputPath;
//#define kOutputPath kStrLocal
//#define kOutputPath kStrFs0Misun
//#define kOutputPath kStrFs0Fuji
//#define kOutputPath kStrTitanJing

//#define Intrepid
//#define Titan
//#define V8
//#define Unknow_machine
extern double CPU_FREQ;
extern char* mach_name;

extern int GROUP_SIZE_IDEAL;
extern int GROUP_SIZE_UPPER_BOUND;

void getfieldname_( int i, char *name );
void getfilename_old(int *id, int *nid);
void adjust_endian();
int swap_int_byte(int *n);
int swap_float_byte(float *n);
int swap_long_long_byte(long long *n);

#endif

