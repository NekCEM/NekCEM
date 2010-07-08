/*  Primary interface for the Phasta Binary read and write routines these*/
/*  functions are 'C' callable.( All arguments have been kept as pointers to*/
/*  facilitate calling from Fortran )*/
/*  Anil Kumar Karanam  Spring 2003*/
#ifdef WIN32
#define openfilerb_ OPENFILERB
#define closefilerb_ CLOSEFILERB
#define readheader_ READHEADER
#define readdatablock_ READDATABLOCK
#define writeheaderrb_ WRITEHEADERRB
#define writedatablock_ WRITEDATABLOCK
#define writestring_ WRITESTRING
#define togglestrictmode_ TOGGLESTRICTMODE
#define SwapArrayByteOrder_ SWAPARRAYBYTEORDER
#define isLittleEndian_ ISLITTLEENDIAN
#ifdef intel
#define for if (0) {} else for
#define ios_base ios
#endif
#elif ( defined IBM )
#define openfilerb_ openfilerb
#define closefilerb_ closefilerb
#define readheader_ readheader
#define readdatablock_ readdatablock
#define writeheaderrb_ writeheaderrb
#define writedatablock_ writedatablock
#define writestring_ writestring
#define togglestrictmode_ togglestrictmode
#endif

#if defined ( __cplusplus )
extern "C" {
#endif
/*added by Jing Fu*/
	void InitPHMPIIO_(int nfields, int nppf, int nfiles, int *filehandle );
	void writerreceive_(int* fileDescriptor, const char datatype[], const char field_name[]);
	void FinalizePHMPIIO_( int *fileDescriptor );
	void setExtra_(int* fileDescriptor, int* writerPos, int* writeMode, const char sharedfileName[]);
/*end of added by Jing Fu*/
    void 
    SwapArrayByteOrder_( void* array, 
                         int   nbytes, 
                         int   nItems ) ;
    void 
    openfilerb_( const char filename[],
               const char mode[],
               int* fileDescriptor );

    void 
    closefile_( int* fileDescriptor,
                const char mode[] );

    void 
    readheader_( int*   fileDescriptor,
                 const char  keyphrase[],
                 void*  valueArray,
                 int*   nItems,
                 const char   datatype[], 
                 const char   iotype[] );
    
    void 
    writeheaderrb_( int*  fileDescriptor,
                  const char keyphrase[],
                  void* valueArray,
                  int*  nItems,
                  int*  ndataItems,
                  const char  datatype[],
                  const char  iotype[] );

    void
    readdatablock_( int*  fileDescriptor,
                    const char keyphrase[],
                    void* valueArray,
                    int*  nItems,
                    const char  datatype[], 
                    const char  iotype[] );
    
    
    void 
    writedatablock_( int*   fileDescriptor,
                     const char  keyphrase[], 
                     void*  valueArray, 
                     int*   nItems, 
                     const char   datatype[], 
                     const char   iotype[]  );


    void 
    writestring_( int* fileDescriptor,
                  const char string[] );

    void
    togglestrictmode_( );

	int
	isLittleEndian_( ) ;

    
#if defined ( __cplusplus )
}
#endif

