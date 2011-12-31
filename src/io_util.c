/**
 * io_util.c contails common functions that is needed by all I/O schemas, e.g.
 * detect endian, swap int, double, float, long long values etc.
 *
 * This file is not dependent on MPI, i.e. they work in both MPI and POSIX
 * environment.
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "io_util.h"

char filename[kMaxPathLen];
int Little_endian = -1;

int DEBUG_FLAG = 0;
int IOTIMER_FLAG = 1;
int IOTRACE_FLAG = 1;

#if defined(Intrepid)
double CPU_FREQ = 850.0 * 1000000;
char* mach_name = "Intrepid@ALCF,ANL";
int GROUP_SIZE_IDEAL = 64;
int GROUP_SIZE_UPPER_BOUND = 64;
#elif defined(Titan)
double CPU_FREQ = 2200.0 * 1000000;
char* mach_name = "Titan@OLCF,ORNL";
int GROUP_SIZE_IDEAL = 64;
int GROUP_SIZE_UPPER_BOUND = 64;
#elif defined(V8)
double CPU_FREQ = 3000.0 * 1000000;
char* mach_name = "V8@MCS,ANL";
int GROUP_SIZE_IDEAL = 64;
int GROUP_SIZE_UPPER_BOUND = 64;
#else
int CPU_FREQ = 2000.0 * 1000000;
char* mach_name = "Unknown machine name (Use 2.0GHz)";
int GROUP_SIZE_IDEAL = 64;
int GROUP_SIZE_UPPER_BOUND = 64;
#endif

/**
 * This function detects machine's endianess and set the global value
 * It's okay if it's called multiple times, should be called at least once
 */
void adjust_endian()
{
	int endian_int = 1;
	char* pchar = (char*) &endian_int;
	if(* (pchar+3)  == 1) Little_endian = 0;
	else Little_endian = 1;
	if(DEBUG_FLAG == 3) {
		if(Little_endian == 0)printf("I'm big endian\n");
		else if(Little_endian == 1) printf("I'm little endian\n");
	}
}

/**
 * This function swap int bytes if it's little endian
 */
int swap_int_byte(int *n)
{
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
	}
	return 0;
}

/**
 * This function swap float bytes if it's little endian
 */
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

/**
 * This function swap double bytes if it's little endian
 * added jingfu 2011-6-29
 */
int swap_double_byte(double *n)
{
	if(Little_endian == 1)
	{
		//printf("entered little endian double swap\n");
		unsigned char *cptr, tmp;
		cptr = (unsigned char*)n;
		tmp = cptr[0];
		cptr[0] = cptr[7];
		cptr[7] = tmp;

		tmp = cptr[1];
		cptr[1] = cptr[6];
		cptr[6] = tmp;

		tmp = cptr[2];
		cptr[2] = cptr[5];
		cptr[5] = tmp;

		tmp = cptr[3];
		cptr[3] = cptr[4];
		cptr[4] = tmp;
	}
	return 0;
}

/**
 * This function swap long long bytes if it's little endian
 * added by Jing Fu at 2010-7-22
 */
int swap_long_long_byte(long long *n)
{
	if(Little_endian == 1)
	{
		unsigned char *cptr, tmp;
		cptr = (unsigned char*)n;
		tmp = cptr[0];
		cptr[0] = cptr[7];
		cptr[7] = tmp;

		tmp = cptr[1];
		cptr[1] = cptr[6];
		cptr[6] = tmp;

		tmp = cptr[2];
		cptr[2] = cptr[5];
		cptr[5] = tmp;

		tmp = cptr[3];
		cptr[3] = cptr[4];
		cptr[4] = tmp;
	}
	return 0;
}

/**
 * This function set field name according to the integer value i
 */
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

/**
 * This function is only to get file name for io_option = 3, i.e. 1PFPP
 */
void getfilename_old(int *id, int *nid )
{
	char ext0[100];
	char ext1[100];

	/*printf( "\n  nid:: %d\n", *nid);*/
	strcpy( filename, "./vtk/binary-NN-p");
	sprintf( ext0, "%.6d-t", *nid);
	strcat( filename, ext0);
	sprintf( ext1, "%.5d", *id);
	strcat( filename, ext1);
	strcat( filename, ".vtk");
	adjust_endian();
}
