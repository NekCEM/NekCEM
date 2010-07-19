#ifndef H_VTKBIN
#define H_VTKBIN

#define ONE_MILLION 1048576

extern char filename[100];
extern char mFilename[100];
extern char rbFilename[100];
extern char rbasciiFilename[100];

extern char* mfBuffer;
extern long long mfBufferCur ;
extern long long mfileCur ;
extern long long fieldSizeSum ;

extern int myrank, mysize;
extern int DEBUG_FLAG;

void getfieldname_( int i, char *name );
void getfilename_(int *id, int *nid );
void adjust_endian();
int swap_int_byte(int *n);
int swap_float_byte(float *n);



#endif

