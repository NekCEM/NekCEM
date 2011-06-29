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

#define cem_out_fields FORTRAN_NAME (cem_out_fields , CEM_OUT_FIELDS )
#define cem_out_fields2 FORTRAN_NAME(cem_out_fields2, CEM_OUT_FIELDS2)
#define cem_out_fields3 FORTRAN_NAME(cem_out_fields3, CEM_OUT_FIELDS3)
#define cem_out_fields4 FORTRAN_NAME(cem_out_fields4, CEM_OUT_FIELDS4)
#define cem_out_fields6 FORTRAN_NAME(cem_out_fields6, CEM_OUT_FIELDS6)

#ifdef UPCASE
void CHECKPOINT_WRITE (int* ioop) 
#elif IBM
void checkpoint_write (int* ioop) 
#else
void checkpoint_write_(int* ioop) 
#endif
{
	if (*ioop == 0) {
		cem_out_fields4();
		cem_out_fields3();
		cem_out_fields2();
	}
	else if (*ioop == 1) {
		cem_out_fields();
	}
	else if (*ioop == 2) {
		cem_out_fields2();
	}
	else if (*ioop == 3) {
		cem_out_fields3();
	}
	else {
#ifdef MPI
		if (*ioop == 4) {
			cem_out_fields4();
		}
		else if (*ioop == 5) {
			set_ascii_nm();
			cem_out_fields6();
		}
		else if (*ioop == 6) {
			cem_out_fields6();
		}
		else if (*ioop == -6) {
			set_ascii_true();
			cem_out_fields6();
		}
		else if (*ioop == 8) {
			set_ascii_nmm();
			cem_out_fields6();
		}
		else {
			printf("ERROR: unknown io_option %d\n", *ioop);
			exit(1);
		}
#endif
	}
}

