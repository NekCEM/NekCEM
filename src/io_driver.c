#ifdef MPI
#include "mpiio_util.h"
#else
#include "io_util.h"
#endif

void checkpoint_write_(int* ioop) {
	if (*ioop == 0) {
		cem_out_fields4_();
		cem_out_fields3_();
		cem_out_fields2_();
	}
	else if (*ioop == 1) {
		cem_out_fields_();
	}
	else if (*ioop == 2) {
		cem_out_fields2_();
	}
	else if (*ioop == 3) {
		cem_out_fields3_();
	}
	else {
#ifdef MPI
		if (*ioop == 4) {
			cem_out_fields4_();
		}
		else if (*ioop == 5) {
			set_ascii_nm_();
			cem_out_fields6_();
		}
		else if (*ioop == 6) {
			cem_out_fields6_();
		}
		else if (*ioop == -6) {
			set_ascii_true_();
			cem_out_fields6_();
		}
		else if (*ioop == 8) {
			set_ascii_nmm_();
			cem_out_fields6_();
		}
		else {
			printf("ERROR: unknown io_option %d\n", *ioop);
			exit(1);
		}
#endif
	}
}

