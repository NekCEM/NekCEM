/* C modules for I/O etc. */
#include <stdio.h>

#ifdef PERFMON
/* #include <perfmon.h> */
#endif

#if   defined NXSRC
#ifndef DELTA
#include <nx.h>
#endif

#ifdef DELTA
#include <mesh.h>
#include <cveclib.h> 
#endif

#elif defined MPISRC
#include <mpi.h>

#endif


extern int my_id;

/*************************************xxt.c************************************
Function: 

Input : 
Output: 
Return: 
Description:  
**************************************xxt.c***********************************/
void
#ifdef UPCASE
HMT_OPCOUNT_START()
#elif  IBM
hmt_opcount_start()
#else
hmt_opcount_start_()
#endif
{
#ifdef PERFMON
  beginmflops();
#endif
}



/*************************************xxt.c************************************
Function: 

Input : 
Output: 
Return: 
Description:  
**************************************xxt.c***********************************/
void
#ifdef UPCASE
HMT_OPCOUNT_END()
#elif  IBM
hmt_opcount_end()
#else
hmt_opcount_end_()
#endif
{
  double opc;

#ifdef PERFMON
  opc = endmflops();
  printmflops();
#endif
}



#ifdef NOT
#include <stdio.h>
#include <perfmon.h>

main()
{
  int i;
  double a[200];
  double sum;
  int rc;
  double drc;

  rc = beginmflops();
  if ( rc != 0 ) {
    perfmonerror( "error calling beginmflops()" );
    exit(1);
  }

  for ( i = 0; i < 200; i++ ) {
      a[i] = i;
  }

  drc = printmflops();
  if ( drc == -1.0 ) {
    perfmonerror( "error calling printmflops()" );
    exit(1);
  }

  for ( i = 0; i < 200; i++ ) {
      sum += a[i];
  }

  drc = endmflops();
  if ( drc == -1.0 ) {
    perfmonerror( "error calling endmflops()" );
    exit(1);
  }

  drc = printmflops();
  if ( drc == -1.0 ) {
    perfmonerror( "error calling printmflops()" );
    exit(1);
  }
}

#endif
