#ifndef H_MPIIO_UTIL
#define H_MPIIO_UTIL
#include <dirent.h>
#include <unistd.h>
#include <sys/stat.h>
#include <pthread.h>

#include <mpi.h>
#include "rdtsc.h"
#include "io_util.h"

#define BGP_FREQ (850*1000000)
#define V8_FREQ (3000*1000000)
#define WRITERBUFFERSIZE (50 * ONE_MILLION)

#define kStrLocal "./vtk"
#define kStrFs0Misun "/intrepid-fs0/users/mmin/scratch/NEKCEM_vtk"
#define kStrFs0Fuji "/intrepid-fs1/users/fuji/scratch/NEKCEM_vtk"

#define kOutputPath kStrLocal
//#define kOutputPath kStrLocal
//#define kOutputPath kStrFs0Misun

extern char rstFilename[100];
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

extern long long start_time;
extern long long end_time;
extern double overall_time;

void getfilename_(int *id, int *nid, int io_option);

#ifdef UPCASE
void SET_ASCII_TRUE ();
#elif  IBM
void set_ascii_true ();
#else
void set_ascii_true_();
#endif

#ifdef UPCASE
void SET_ASCII_NMM ();
#elif  IBM
void set_ascii_nmm ();
#else
void set_ascii_nmm_();
#endif

#ifdef UPCASE
void SET_ASCII_NM ();
#elif  IBM
void set_ascii_nm ();
#else
void set_ascii_nm_();
#endif

// define file structure to pass to IO thread
typedef struct file_s {
	pthread_t pthread;
	pthread_attr_t thread_attr;
	pthread_mutex_t mutex;
	pthread_mutexattr_t mutex_attr;
	pthread_cond_t	cond;

	MPI_File* pmfile; // pointer to file descriptor
//	long long* pmfileCur; // pointer to file current position
	char* pwriterBuffer; // pointer to write buffer
	long long llwriterBufferCur; // pointer to write buffer current position
} file_t;
// following for rbIO_nekcem.c
void init_file_struc();
void reset_file_struc();
void free_file_struc(file_t* file);
void run_io_thread(file_t* file);

#endif

