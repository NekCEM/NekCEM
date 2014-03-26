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

//###############################################################################
// if prod_mode is defined, there will be no debug info print out;
// otherwise (if commented out), debug will be on, and i/o will happen in local dir
//###############################################################################
#define prod_mode
//###############################################################################

char filename[kMaxPathLen];
int Little_endian = -1;

//int DEBUG_FLAG = 1;
int IOTIMER_FLAG = 1;
int IOTRACE_FLAG = 0;       // off by default, set by param(86) in *.rea file TODO: pass param from cem_dg2 to I/O, disabled
int COMPUTE_TRACE_FLAG = 0; // off by default, set by param(85) in *.rea file TODO: pass param from cem_dg2 to I/O, disabled

int io_step; // global io_step number

// these def's are in makenek script so to be automatically machine specific
// you have to think really hard before changing anything here..
#if defined(Intrepid)
double CPU_FREQ = 850.0 * 1000000;
int GROUP_SIZE_IDEAL = 64;
int GROUP_SIZE_UPPER_BOUND = 64;
long long WRITERBUFFERSIZE = (200 * ONE_MILLION);
char* mach_name = "Intrepid@ALCF,ANL";
#ifdef prod_mode
//char* kOutputPath = "vtk"; //please note, as of 03/25/2012, this sets hierarchical setting for output files with subdirectories under local ./vtk
char* kOutputPath = kStrLocal; //please note, as of 03/15/2012, intrepid home fs is read-only; so compile/run code on fs0 if you want to use local
//char* kOutputPath = kStrFs0Fuji;  // This is used for Jing's test purpose
int DEBUG_FLAG = 0; // no io-print statement
#else
char* kOutputPath = kStrLocal;
int DEBUG_FLAG = 1;
#endif

#elif defined(Titan)
double CPU_FREQ = 2200.0 * 1000000;
int GROUP_SIZE_IDEAL = 64;
int GROUP_SIZE_UPPER_BOUND = 64;
long long WRITERBUFFERSIZE = (200 * ONE_MILLION);
char* mach_name = "Titan@OLCF,ORNL";
#ifdef prod_mode
//char* kOutputPath = "vtk"; ////please note, as of 03/25/2012, this sets hierarchical setting for output files with subdirectories under local ./vtk
char* kOutputPath = kStrLocal;
//char* kOutputPath = kStrTitanJing;
int DEBUG_FLAG = 0; // no io-print statement
#else
char* kOutputPath = kStrLocal;
int DEBUG_FLAG = 1;
int DEBUG_FLAG = 0; // no io-print statement
#endif

#else
double CPU_FREQ = 2000.0 * 1000000;
int GROUP_SIZE_IDEAL = 64;
int GROUP_SIZE_UPPER_BOUND = 64;
long long WRITERBUFFERSIZE = (50 * ONE_MILLION);
char* kOutputPath = kStrLocal;
char* mach_name = "Unknown machine name (Use 2.0GHz)";

int DEBUG_FLAG = 0;
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
                        strcpy(name, "E ");
                        break;
                case 2:
                        strcpy(name, "H ");
                        break;
                case 3:
                        strcpy(name, "solE");
                        break;
                case 4:
                        strcpy(name, "solH");
                        break;
                case 11:
                        strcpy(name, "sch_Ur_Ui_Ui");
                        break;
                case 12:
                        strcpy(name, "sch_Qr");
                        break;
                case 13:
                        strcpy(name, "sch_Qi");
                        break;
                case 14:
                        strcpy(name, "sch_sUr_sQr_sQi");
                        break;
                case 21:
                        strcpy(name, "JN");
                        break;
                case 22:
                        strcpy(name, "PN");
                        break;
                case 23:
                        strcpy(name, "QN");
                        break;
                case 24:
                        strcpy(name, "sJN");
                        break;
                case 31:
                        strcpy(name, "cN_cP_cE");
                        break;
                case 32:
                        strcpy(name, "scN_scP_scE");
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
