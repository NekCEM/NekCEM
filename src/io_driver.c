/**
 * This file contains the high-level functions that decides which checkpoint schema to use
 * according to the io_option integer value.
 *
 * This can be used in either MPI or non-MPI environment.
 */

#ifdef MPI
#include "mpiio_util.h"
#else
#include "io_util.h"
#endif
//#include "jl2/name.h"

#ifdef UPCASE
#  define FORTRAN_NAME(low,up) up
#else
#ifdef UNDERSCORE
#  define FORTRAN_NAME(low,up) low##_
#else
#  define FORTRAN_NAME(low,up) low
#endif
#endif

#define cem_restart_out FORTRAN_NAME(cem_restart_out, CEM_RESTART_OUT)
#define cem_out_fields  FORTRAN_NAME(cem_out_fields , CEM_OUT_FIELDS )
#define cem_out_fields2 FORTRAN_NAME(cem_out_fields2, CEM_OUT_FIELDS2)
#define cem_out_fields3 FORTRAN_NAME(cem_out_fields3, CEM_OUT_FIELDS3)
#define cem_out_fields4 FORTRAN_NAME(cem_out_fields4, CEM_OUT_FIELDS4)
#define cem_out_fields6 FORTRAN_NAME(cem_out_fields6, CEM_OUT_FIELDS6)

#ifdef UPCASE
void CHECKPOINT_WRITE (int* ioop, int* idump)
#elif IBM
void checkpoint_write (int* ioop, int* idump)
#else
void checkpoint_write_(int* ioop, int* idump)
#endif
{

  //if(DEBUG_FLAG) printf("in checkpoint_write()\n");
  io_step ++;
	if       (*ioop == 0) {
		cem_out_fields4(idump);
		cem_out_fields3(idump);
		cem_out_fields2(idump);
	}
	else if (*ioop == 1) {
		cem_out_fields (idump);
	}
	else if (*ioop == 2) {
		cem_out_fields2(idump);
	}
	else if (*ioop == 3) {
		cem_out_fields3(idump);
	}
	else {
#ifdef MPI
		if      (*ioop == 99) {
			cem_restart_out(idump);
		}
		else if (*ioop == 4) {
			cem_out_fields4(idump);
		}
		else if (*ioop == 5) {
      set_io_option(5);
			//set_ascii_nm();
			cem_out_fields6(idump);
		}
		else if (*ioop == 6) {
      set_io_option(6);
			cem_out_fields6(idump);
		}
		else if (*ioop ==-6) {
      set_io_option(-6);
			//set_ascii_true();
			cem_out_fields6(idump);
		}
		else if (*ioop == 8) {
      set_io_option(8);
			//set_ascii_nmm();
			cem_out_fields6(idump);
		}
		else if (*ioop == 18) {
      set_io_option(18);
			//set_ascii_nmm();
			cem_out_fields6(idump);
		}
		else {
			printf("ERROR: unknown io_option %d\n", *ioop);
			exit(1);
		}
#endif
	}
}

