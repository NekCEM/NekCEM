#ifndef H_VTKBIN
#define H_VTKBIN

void getfieldname_( int i, char *name );
void getfilename_(int *id, int *nid );
void adjust_endian();
int swap_int_byte(int *n);
int swap_float_byte(float *n);



#endif

