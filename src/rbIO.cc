#include <map>
#include <vector>
#include <string>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>
#include <new>
#include <iostream>
#include "rbIO.h"
#include "mpi.h"
//#include "rdtsc.h"

#define VERSION_INFO_HEADER_SIZE 2048
#define DB_HEADER_SIZE   1024
#define ONE_MEGABYTE  1048576
#define ENDIAN_TEST_NUMBER 12180 // Troy's Zip Code!!
#define MAX_PHASTA_FILES 1024
#define MAX_PHASTA_FILE_NAME_LENGTH 1024
#define MAX_FIELDS_NUMBER 16
#define MAX_FIELDS_NAME_LENGTH 128
#define MASTER_HEADER_SIZE ONE_MEGABYTE

//####################
//for rbIO
//####################
#define MAX_WRITER_BUFFER_SIZE_MB 250*ONE_MEGABYTE

//end of for rbIO

enum PhastaIO_Errors
  {
    MAX_PHASTA_FILES_EXCEEDED = -1,
    UNABLE_TO_OPEN_FILE = -2,
    NOT_A_MPI_FILE = -3,
    GPID_EXCEEDED = -4,
    DATA_TYPE_ILLEGAL = -5,
  };

using namespace std;
//using namespace rbIO;
#define swap_char(A,B) { ucTmp = A; A = B ; B = ucTmp; }

// now we create a local "unnamed" namespace which will contain all the
// functions which will be used locally within this file.  These functions
// donot have external declerations and are not visible outside this file.

int rdtsc(){return 200;}

namespace{
    
  map< int , char* > LastHeaderKey;
  vector< FILE* > fileArray;
  vector< bool > byte_order;
  vector< int > header_type;
  int DataSize=0;
  bool LastHeaderNotFound = false;
  bool Wrong_Endian = false ;
  bool Strict_Error = false ;
  bool binary_format = true;

  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //not useful, just for compiling
 
  vector< MPI_File > MPI_Handle;
  bool Parallel_IO = true;
  int file_header_offset = ONE_MEGABYTE;  //the place where you start writing db_header
  int datablock_offset = DB_HEADER_SIZE;  // this is the size of a db_header
  int MPI_Tag_space = VERSION_INFO_HEADER_SIZE;   // this is the size of version_info_header_size

  //#################################
  // rbIO code starts from here
  //################################


  typedef struct{
    char key[128];
    unsigned long long offset;
    int GPid;
  } header_array_t;

  // bool Wrong_Endian;                            /* default to false */
  typedef struct
  {
    bool Wrong_Endian;
    
    char fileName[MAX_PHASTA_FILE_NAME_LENGTH];   /* defafults to 1024 */
    char sharedfileName[MAX_PHASTA_FILE_NAME_LENGTH];//for shared file, mode 1 
    int nppp;                                   
    int nPPF;                                   
    int nFiles;                                 
    int nFields;                                
    long long my_offset;//looks like it's never used in rbIO writing
    MPI_File file_handle; //not used

    double* master_header;

    double * double_chunk[MAX_FIELDS_NUMBER];
    int * int_chunk[MAX_FIELDS_NUMBER];

    double * read_double_chunk;
    int * read_int_chunk;

    long long **my_offset_table;
    unsigned long long **my_read_table;

    int field_count ;
    int part_count;
    int read_field_count; 
    int read_part_count;
    int GPid; 
    int start_id;
    unsigned long long next_start_address;         /* default to file_header_offset */

    //header_array_t **offset_array;
    //split the comm_world here:
    int myrank;
    int numprocs;
    int local_myrank;
    int local_numprocs;
    MPI_Comm local_comm;

    //for rbIO

    char field_table[MAX_FIELDS_NUMBER][MAX_FIELDS_NAME_LENGTH];//hold all fields names

    long long totalFileSize ;
    long long thisFieldSize ;
    long long data_location;
 
    double* writer_double_buffer; //the big buffer to hold all data
    int * writer_int_buffer;
    int writer_datatype; // "1" stands for "double", "2" for "integer"
  
    MPI_File mfile;
    char datablock_header[DB_HEADER_SIZE];	
   
    int commSize, rankInComm, numGroups, groupSize, groupRank, rankInGroup;
    int writerPos ;
    int realwriterPos ;  //after writerpos transform, realwriterpos is here for rank computing
      
    int writeMode; //default to 0, meaning every writer writes its own file; 1 means to shared file
	
 
    //debug memory management problems
    //double * temp_first;
    //double * allData;
    //double * field_data_send;

    //timing variables
    long long *irecv_start, *irecv_end;
    double *irecv_time;
    
  } rbio_file_t;

  rbio_file_t *rbIOFileA[MAX_PHASTA_FILES];
  int rbIONextIndex = 0; //indicates next index to allocate

  bool DEBUG_FLAG, TIMING_FLAG, AUGMENT_FLAG;


  char*
  StringStripper( const char  istring[] ) {

    int length = strlen( istring );
    char* dest = new char [ length + 1 ];
    strcpy( dest, istring );
    dest[ length ] = '\0';

    if ( char* p = strpbrk( dest, " ") ) 
      *p = '\0';
	
    return dest;
  }

  inline int 
  cscompare( const char teststring[], 
	     const char targetstring[] ) {

    char* s1 = const_cast<char*>(teststring);
    char* s2 = const_cast<char*>(targetstring);

    while( *s1 == ' ') s1++;
    while( *s2 == ' ') s2++;
    while( ( *s1 ) 
	   && ( *s2 ) 
	   && ( *s2 != '?')
	   && ( tolower( *s1 )==tolower( *s2 ) ) ) {
      s1++;
      s2++;
      while( *s1 == ' ') s1++;
      while( *s2 == ' ') s2++;
    }
    if ( !( *s1 ) || ( *s1 == '?') ) return 1;
    else return 0;
  }
    
  inline void
  isBinary( const char iotype[] ) {

    char* fname = StringStripper( iotype );
    if ( cscompare( fname, "binary" ) ) binary_format = true;
    else binary_format = false;
    delete [] fname;

  }

  inline size_t
  typeSize( const char typestring[] ) {

    char* ts1 = StringStripper( typestring );

    if ( cscompare( "integer", ts1 ) ) {
      delete [] ts1;
      return sizeof(int);
    } else if ( cscompare( "double", ts1 ) ) { 
      delete [] ts1;
      return sizeof( double );
    } else { 
      delete [] ts1;
      fprintf(stderr,"unknown type : %s\n",ts1);
      return 0;
    }
  }
    
  int
  readHeader( FILE*       fileObject,
	      const char  phrase[],
	      int*        params,
	      int         expect ) {
        
    char* text_header;
    char* token;
    char Line[1024];
    char junk;
    bool FOUND = false ;
    int real_length;
    unsigned int skip_size;
    int integer_value;
    int rewind_count=0;

    if( !fgets( Line, 1024, fileObject ) && feof( fileObject ) ) {
      rewind( fileObject );
      clearerr( fileObject );
      rewind_count++;
      fgets( Line, 1024, fileObject );
    }


    while( !FOUND  && ( rewind_count < 2 ) )  {
      if ( ( Line[0] != '\n' ) && ( real_length = strcspn( Line, "#" )) ){
	text_header = new char [ real_length + 1 ];
	strncpy( text_header, Line, real_length );
	text_header[ real_length ] =static_cast<char>(NULL);
	token = strtok ( text_header, ":" );
	if( cscompare( phrase , token ) ) {
	  FOUND = true ;
	  token = strtok( NULL, " ,;<>" );	     
	  skip_size = atoi( token );
	  int i;
	  for( i=0; i < expect && ( token = strtok( NULL," ,;<>") ); i++) {
	    params[i] = atoi( token );
	  }
	  if ( i < expect ) {
	    fprintf(stderr,"Expected # of ints not found for: %s\n",phrase );
	  }
	} else if ( cscompare(token,"byteorder magic number") ) {
	  if ( binary_format ) {
	    fread((void*)&integer_value,sizeof(int),1,fileObject);
	    fread( &junk, sizeof(char), 1 , fileObject );
	    if ( 362436 != integer_value ) Wrong_Endian = true;
	  } else{
	    fscanf(fileObject, "%d\n", &integer_value );
	  }
	} else { 
	  /* some other header, so just skip over */

	  token = strtok( NULL, " ,;<>" );
	  skip_size = atoi( token );
	  if ( binary_format) 
	    fseek( fileObject, skip_size, SEEK_CUR );
	  else 
	    for(unsigned int gama=0; gama < skip_size; gama++ ) 
	      fgets( Line, 1024, fileObject );
	}
	delete [] text_header;
      }
	

      if ( !FOUND ) 
	if( !fgets( Line, 1024, fileObject ) && feof( fileObject ) ) {
	  rewind( fileObject );
	  clearerr( fileObject );
	  rewind_count++;
	  fgets( Line, 1024, fileObject );
	}
    }             
        
    if ( !FOUND ) {
      fprintf(stderr, "Error: Cound not find: %s\n", phrase);
      return 1;
    }
    return 0;
  }	

} // end unnamed namespace
 

void HelloWorld()
{
printf("Hello World\n");
}
 
// begin of publicly visible functions
void InitPHMPIIO_( int nfields, int nppf, int nfiles, int *filehandle )
{
 
  if( rbIONextIndex == MAX_PHASTA_FILES)
    {
      //return MAX_PHASTA_FILES_EXCEEDED;
    }
  else if( rbIONextIndex == 0)
    {
     // for( int i = 0 ; i < MAX_PHASTA_FILES; i++)
	//rbIOFileA[i] = NULL;
    }

  
  rbIOFileA[rbIONextIndex] = (rbio_file_t*) calloc(1, sizeof(rbio_file_t) );
  int i = rbIONextIndex ;
  rbIONextIndex ++;
  *filehandle = i;// the library takes charge to increment file descriptor, users do nothing
  
  rbIOFileA[i]->Wrong_Endian = false;


  rbIOFileA[i]-> numGroups = nfiles;
  rbIOFileA[i]->nFiles = nfiles;
  rbIOFileA[i]->nPPF = nppf;
  rbIOFileA[i]->nFields = nfields;	
  rbIOFileA[i]->nppp = 1;

  MPI_Comm_size(MPI_COMM_WORLD, & (rbIOFileA[i]->commSize) );
  MPI_Comm_rank(MPI_COMM_WORLD, & (rbIOFileA[i]->rankInComm) );
//rbIOFileA[i]->rankInComm = workfc.globalrank;
//rbIOFileA[i]->commSize = workfc.globalsize;

	
  rbIOFileA[i]->groupSize =   (rbIOFileA[i]->commSize) / ( rbIOFileA[i]->numGroups );
  rbIOFileA[i]->groupRank = ( rbIOFileA[i]-> rankInComm)/(  rbIOFileA[i]->groupSize);
  rbIOFileA[i]->rankInGroup = (  rbIOFileA[i]->rankInComm) % (  rbIOFileA[i]->groupSize);

  DEBUG_FLAG = false ;
  TIMING_FLAG = true;
  AUGMENT_FLAG = true;

  rbIOFileA[i]->writeMode = 0; //default to 0

  //**********************************
  //writerPos = *filehandle;
  //no way to pass writerPos value to the lib yet, add manually
  rbIOFileA[i]->writerPos = 0;
  //***********************************


  rbIOFileA[i]->realwriterPos = rbIOFileA[i]-> writerPos;

  //do the transform, make non-zero writer pos to zero
  //

  if(rbIOFileA[i]->rankInGroup < rbIOFileA[i]-> writerPos)
    rbIOFileA[i]->rankInGroup ++;
  else if(rbIOFileA[i]->rankInGroup == rbIOFileA[i]->writerPos)
    {
      rbIOFileA[i]->rankInGroup = 0;
    }
  rbIOFileA[i]->writerPos = 0;

  if(rbIOFileA[i]->rankInGroup == rbIOFileA[i]->writerPos)
    {
//     rbIOFileA[i]->my_offset_table = new long long*[rbIOFileA[i]->nFields];
	rbIOFileA[i]->my_offset_table = ( long long**) (sizeof( long long*) * rbIOFileA[i]->nFields );
     for( int j = 0; j <  rbIOFileA[i]->nFields; j++)
       {
//	 rbIOFileA[i]->my_offset_table[j] = new long long[ rbIOFileA[i]->groupSize -1 ];
	rbIOFileA[i]->my_offset_table[j] = (long long*) ( sizeof(long long ) * (rbIOFileA[i]->groupSize -1) ) ;

       }
//DEBUG_FLAG= true;////////////////////////////////////////////////////////////////
  }

/*
  if(DEBUG_FLAG)
    {
      printf("myrank = %d, commSize = %d, my groupRank = %d, my rankInGroup = %d,groupSize = %d, writerPos is %d\n", 
	     rbIOFileA[i]->rankInComm,
	     rbIOFileA[i]->commSize, 
	     rbIOFileA[i]->groupRank, 
	     rbIOFileA[i]->rankInGroup, 
	     rbIOFileA[i]->groupSize, 
	     rbIOFileA[i]->writerPos);
    }
*/

  return ;	
}

void setExtra_(int* fileDescriptor,
	      int* writerPos,
	      int* writeMode,
	      const char filename[])
{
	int i = *fileDescriptor;
	int myColor = -1;

	rbIOFileA[i]->writeMode = *writeMode;
	strcpy(rbIOFileA[i]->sharedfileName, filename);

	//set writer_comm
	if(rbIOFileA[i]->rankInGroup == rbIOFileA[i]->writerPos)
		myColor = 1;
	else	
		myColor = 2;
	
	MPI_Comm_split(MPI_COMM_WORLD, myColor, rbIOFileA[i]->rankInComm, &(rbIOFileA[i]->local_comm));

	
}
void FinalizePHMPIIO_( int *fileDescriptor )
{
    int i, j;
    i = *fileDescriptor;

//added by Jing Fu
rbIONextIndex = 0;
//add end

    if(rbIOFileA[i]->rankInGroup == rbIOFileA[i]->writerPos)
      {
	for( int j = 0; j <  rbIOFileA[i]->nFields; j++)
	  {
	    delete [] rbIOFileA[i]->my_offset_table[j] ;
	  }
	delete [] rbIOFileA[i]->my_offset_table;
      }
    delete [] rbIOFileA[i]->irecv_start;
    delete [] rbIOFileA[i]->irecv_end;
    delete [] rbIOFileA[i]->irecv_time;

    free( rbIOFileA[i]); //this was alloc by "calloc"
}


void 
SwapArrayByteOrder_( void* array, 
                     int   nbytes, 
                     int   nItems ) {
  /* This swaps the byte order for the array of nItems each
     of size nbytes , This will be called only locally  */
  int i,j;
  unsigned char ucTmp;
  unsigned char* ucDst = (unsigned char*)array;
        
  for(i=0; i < nItems; i++) {
    for(j=0; j < (nbytes/2); j++)
      swap_char( ucDst[j] , ucDst[(nbytes - 1) - j] );
    ucDst += nbytes;
  }
}
    
void
openfilerb_( const char filename[],
	   const char mode[],
	   int*  fileDescriptor ) {

  int j = *fileDescriptor;
  
  if( rbIONextIndex == 0)//it's serial lib
    {
      FILE* file=NULL ;
      *fileDescriptor = 0;
      char* fname = StringStripper( filename );
      char* imode = StringStripper( mode );

      if ( cscompare( "read", imode ) ) file = fopen(fname, "rb" );
      else if( cscompare( "write", imode ) ) file = fopen(fname, "wb" );
      else if( cscompare( "append", imode ) ) file = fopen(fname, "ab" );

    
      if ( !file ){

	fprintf(stderr,"unable to open file : %s\n",fname ) ;

      } else {

	fileArray.push_back( file );
	byte_order.push_back( false );         
	header_type.push_back( sizeof(int) );
	*fileDescriptor = fileArray.size();
        
      }
      delete [] fname;
      delete [] imode;
     // printf("serial openfile_ done: %s - %s\n", filename, mode);/////////////////////////////////////////////////////
    }
  else //it's parallel io
    {

      char* fname = StringStripper( filename );
      char* imode = StringStripper( mode );
      strcpy( rbIOFileA[j]->fileName, filename);

      if ( cscompare( "write", imode ) ){
	if(rbIOFileA[j]->rankInGroup == rbIOFileA[j]->writerPos && rbIOFileA[j]->writeMode == 0){
	  //here I'll firstly detect if there exist files with the same name,
	  //if there are, delete first, then re-create 
	  //this is to make sure performance test accurate and also avoid file append problems

	  //printf("writer deleting previous files and create new ones\n");
	  int temp = MPI_File_open(MPI_COMM_SELF, rbIOFileA[j]->fileName, MPI_MODE_CREATE|MPI_MODE_EXCL, MPI_INFO_NULL, &(rbIOFileA[j]->mfile));
	  if( temp)
	    {
	      //                      printf("file exists...deleting..\n");
	      int revdlt = MPI_File_delete(rbIOFileA[j]->fileName, MPI_INFO_NULL);
	      //                      if(revdlt)printf("delete original file failed!");
	    }

	  int rc;

	  rc = MPI_File_open(MPI_COMM_SELF, rbIOFileA[j]->fileName, MPI_MODE_CREATE | MPI_MODE_RDWR , MPI_INFO_NULL, &(rbIOFileA[j]->mfile));
	  if(rc){
	    printf("Unable to create file %s \n", rbIOFileA[j]->fileName);
	    fflush(stdout);
	  }

	  //########################################################
	  //should re-set field_count and free double chunk etc???
	  //or make sure everything is set in close_file function
	  //######################################################

	  rbIOFileA[j]->field_count = 0;
	  rbIOFileA[j]->writer_datatype = -1;
	  rbIOFileA[j]->totalFileSize = MASTER_HEADER_SIZE;

	 printf("writer created new files\n");
	}

	/*******if it's writer and it's shared file**********/
	/****************************************************/
	else if(rbIOFileA[j]->rankInGroup == rbIOFileA[j]->writerPos && rbIOFileA[j]->writeMode == 1){
	//use writer0 to detect and delete old files
	if(rbIOFileA[j]->groupRank == 1)
	{
          int temp = MPI_File_open(MPI_COMM_SELF, rbIOFileA[j]->sharedfileName, MPI_MODE_CREATE|MPI_MODE_EXCL, MPI_INFO_NULL, &(rbIOFileA[j]->mfile));
          if( temp)
            {
              //                      printf("file exists...deleting..\n");
              int revdlt = MPI_File_delete(rbIOFileA[j]->sharedfileName, MPI_INFO_NULL);
              //                      if(revdlt)printf("delete original file failed!");
            }
	}
	MPI_Barrier(rbIOFileA[j]->local_comm);

          int rc;
          rc = MPI_File_open(rbIOFileA[j]->local_comm, rbIOFileA[j]->sharedfileName, MPI_MODE_CREATE | MPI_MODE_RDWR , MPI_INFO_NULL, &(rbIOFileA[j]->mfile));
          if(rc){
            printf("Unable to create shared file %s \n", rbIOFileA[j]->sharedfileName);
            fflush(stdout);
          }

          //########################################################
          //should re-set field_count and free double chunk etc???
          //or make sure everything is set in close_file function
          //######################################################

          rbIOFileA[j]->field_count = 0;
          rbIOFileA[j]->writer_datatype = -1;
          rbIOFileA[j]->totalFileSize = MASTER_HEADER_SIZE;

        }

      //printf(" parallel open file done: %s - %s\n", filename, mode);//////////////////////////////////////////////////////
      }//end of if it's "write"
      else if  ( cscompare( "read", imode ) )
	{
	}
      else
	{
	  //error
	}
    }//end of if parallel io
}

void
closefilerb_( int* fileDescriptor,
	       const char mode[] ) {

  int j = *fileDescriptor;

  if( rbIONextIndex == 0) //it's serial lib
    {
      
      char* imode = StringStripper( mode );

      if( cscompare( "write", imode ) 
	  || cscompare( "append", imode ) ) {
	fflush( fileArray[ *fileDescriptor - 1 ] );
      } 

      fclose( fileArray[ *fileDescriptor - 1 ] );
      delete [] imode;
//      printf("serial close done: file-%d\n", *fileDescriptor) ;/////////////////////////////////////////////////
    }
  else //parallel lib
    {
      if( cscompare( "write", mode))
	{
	  //if it is writer, need to flush master header and offset table before closing file
	  if( rbIOFileA[j]->rankInGroup ==  rbIOFileA[j]->writerPos)
	    {
	      //########################################
	      //write version_info
	      //#############################################

	      //rbIOFileA[j]->master_header = new double[MASTER_HEADER_SIZE/sizeof(double)];
	
	      char master_header_version_info[VERSION_INFO_HEADER_SIZE];
	      bzero((void*)master_header_version_info, VERSION_INFO_HEADER_SIZE);

	      sprintf(master_header_version_info, "MPI_IO_TAG :");
	  
	      //have this bizzarre thing coz' Ning is using seperate integer read/write to magic number	 
	      int write_magic_number = ENDIAN_TEST_NUMBER;
	      memcpy(&master_header_version_info[sizeof("MPI_IO_TAG :")], &write_magic_number, sizeof(int));

	      sprintf(&master_header_version_info[MAX_FIELDS_NAME_LENGTH], "\nnFields : %d\n",  rbIOFileA[j]->nFields);
	      for( int i = 0; i <  rbIOFileA[j]->field_count; i++)
		{
		  sprintf(&master_header_version_info[ (2+i)*MAX_FIELDS_NAME_LENGTH ], "\n%s", & rbIOFileA[j]->field_table[i]); 
		}

	      sprintf(&master_header_version_info[(2+ rbIOFileA[j]->field_count)*MAX_FIELDS_NAME_LENGTH], "\nnPPF : %d\n",  rbIOFileA[j]->nPPF);
	      //sprintf(&master_header_version_info[(3+ rbIOFileA[j]->field_count)*MAX_FIELDS_NAME_LENGTH], "nFiles:%d",  rbIOFileA[j]->nFiles);

	      //add size info for reading convenience
	      if(AUGMENT_FLAG)
	      memcpy(&master_header_version_info[VERSION_INFO_HEADER_SIZE-sizeof(long long)], & rbIOFileA[j]->totalFileSize, sizeof(long long));

	      //for(int k = 0; k < VERSION_INFO_HEADER_SIZE; k++)
	      //printf("%c-", master_header_version_info[k]);

	      //#######################################
	      //if data type is double for writer
	      //######################################
	      if( rbIOFileA[j]->writer_datatype == 1)
		{
	      
		  memcpy( rbIOFileA[j]->writer_double_buffer, master_header_version_info, VERSION_INFO_HEADER_SIZE);

		  //write offset_table
		  //printf("my field_count is %d\n",  rbIOFileA[j]->field_count);
		  for(int i = 0; i <  rbIOFileA[j]->field_count; i++)
		    {
		      memcpy(& (rbIOFileA[j]->writer_double_buffer[(VERSION_INFO_HEADER_SIZE + i*( rbIOFileA[j]->groupSize-1)*sizeof(long long)) / sizeof(double)]), 
			     rbIOFileA[j]->my_offset_table[i],
			     ( rbIOFileA[j]->groupSize-1)*sizeof(long long));


		     
		      
		      if(DEBUG_FLAG)
			{
			  for( int k = 0; k <  rbIOFileA[j]->groupSize -1; k++) rbIOFileA[j]->my_offset_table[i][k] = 0;
			  memcpy( rbIOFileA[j]->my_offset_table[ i ], 
				  & (rbIOFileA[j]->writer_double_buffer[(VERSION_INFO_HEADER_SIZE + i*( rbIOFileA[j]->groupSize-1)*sizeof(long long)) / sizeof(double)]), 
				  ( rbIOFileA[j]->groupSize-1)*sizeof(long long));
			  for( int k = 0; k <  rbIOFileA[j]->groupSize -1; k++) 
			    {
	   
			      printf("++ %lld = \n",  rbIOFileA[j]->my_offset_table[ i ][k]);
			    }
			}
	
		    }
	  
		  if(DEBUG_FLAG)
			printf("master header done..starting memcpy and write\n");
		  if(DEBUG_FLAG ) 
			printf("before File_write_at(), totalFileSize is %lld\n",  rbIOFileA[j]->totalFileSize);
	

		  MPI_Status write_data_status;
		  if(rbIOFileA[j]->writeMode == 0)
		{
		  MPI_File_write_at( rbIOFileA[j]->mfile,
				     0, 
				     & (rbIOFileA[j]->writer_double_buffer[0]),  
				     rbIOFileA[j]->totalFileSize/sizeof(double), 
				     MPI_DOUBLE, 
				     &write_data_status);

		  MPI_File_close( & (rbIOFileA[j]->mfile) );
	      	}
		  else if( rbIOFileA[j]->writeMode == 1)
		{
		  MPI_Barrier(rbIOFileA[j]->local_comm);
		  long long write_offset;
		  MPI_Scan(&rbIOFileA[j]->totalFileSize, &write_offset, 1, MPI_LONG_LONG, MPI_SUM, rbIOFileA[j]->local_comm);
		  MPI_File_write_at_all_begin( rbIOFileA[j]->mfile,
						write_offset,
						& (rbIOFileA[j]->writer_double_buffer[0]),
                                     		rbIOFileA[j]->totalFileSize/sizeof(double),
                                     		MPI_DOUBLE);
		  MPI_File_write_at_all_end( 	rbIOFileA[j]->mfile,
						& (rbIOFileA[j]->writer_double_buffer[0]),
                                     		&write_data_status);
		  MPI_File_close( & (rbIOFileA[j]->mfile) );

		}
		  for(int i = 0; i < MAX_FIELDS_NUMBER ; i++)
		    {
		      delete [] rbIOFileA[j]->double_chunk[i];
		    }
		  delete [] rbIOFileA[j]->writer_double_buffer;
		  // delete [] rbIOFileA[j]->master_header;
	      
		}
	      //---------------------------------------------------------------------------------
	      //#######################################
	      //if data type is int for writer
	      //######################################
	      else if(rbIOFileA[j]->writer_datatype == 2)//it's integer
		{
		  memcpy( rbIOFileA[j]->writer_int_buffer, master_header_version_info, VERSION_INFO_HEADER_SIZE);

		  //write offset_table
		  for(int i = 0; i <  rbIOFileA[j]->field_count; i++)
		    {
		      memcpy(& (rbIOFileA[j]->writer_int_buffer[(VERSION_INFO_HEADER_SIZE + i*( rbIOFileA[j]->groupSize-1)*sizeof(long long))/sizeof(int)]), 
			     rbIOFileA[j]->my_offset_table[i],
			     ( rbIOFileA[j]->groupSize-1)*sizeof(long long));
		    }
	  
		  if(DEBUG_FLAG)printf("master header done..starting memcpy and write\n");
	      
		  if(DEBUG_FLAG ) printf("before File_write_at(), totalFileSize is %lld\n",  rbIOFileA[j]->totalFileSize);

		  MPI_Status write_data_status;
		  MPI_File_write_at( rbIOFileA[j]->mfile,
				     0, 
				     & (rbIOFileA[j]->writer_int_buffer[0]),  
				     rbIOFileA[j]->totalFileSize/sizeof(int), 
				     MPI_INT, 
				     &write_data_status);

		  MPI_File_close( & (rbIOFileA[j]->mfile) );
	      
		  for(int i = 0; i < MAX_FIELDS_NUMBER ; i++)
		    {
		      delete [] rbIOFileA[j]->int_chunk[i];
		    }
		  delete [] rbIOFileA[j]->writer_int_buffer;
		}
	    }//end of if writer
	  //---------------------------------------------------------------------------------
	  else  //it's a worker, now it's ok to free the send buffer
	    {
	      if(rbIOFileA[j]->writer_datatype == 1)
		{
		  for(int i = 0; i <  rbIOFileA[j]->field_count ; i++)
		    {
		      delete [] (rbIOFileA[j]->double_chunk[i]);
		    }
		}
	      else if(rbIOFileA[j]->writer_datatype == 2)
		{
		  for(int i = 0; i <  rbIOFileA[j]->field_count ; i++)
		    {
		      delete [] (rbIOFileA[j]->int_chunk[i]);
		    }
		}
	      
	    }//end of if worker
	//printf("parallel closefile_ done\n");
	}//end of if "write"
    }
}


void
readheader_( int* fileDescriptor,
             const char keyphrase[],
             void* valueArray,
             int*  nItems,
             const char  datatype[],
             const char  iotype[] ) {

  int i = *fileDescriptor;

  if(rbIONextIndex == 0)//serial lib
    { 
      int filePtr = *fileDescriptor - 1;
      FILE* fileObject;
      int* valueListInt;

      if ( *fileDescriptor < 1 || *fileDescriptor > (int)fileArray.size() ) {
	fprintf(stderr,"No file associated with Descriptor %d\n",*fileDescriptor);
	fprintf(stderr,"openfile_ function has to be called before \n") ;
	fprintf(stderr,"acessing the file\n ") ;
	fprintf(stderr,"fatal error: cannot continue, returning out of call\n");
	return;
      }

      LastHeaderKey[ filePtr ] = const_cast< char* >( keyphrase ); 
      LastHeaderNotFound = false;

      fileObject = fileArray[ filePtr ] ;
      Wrong_Endian = byte_order[ filePtr ];

      isBinary( iotype );
      typeSize( datatype );   //redundant call, just avoid a compiler warning.

      // right now we are making the assumption that we will only write integers
      // on the header line.

      valueListInt = static_cast< int* >( valueArray );
      int ierr = readHeader( fileObject ,
			     keyphrase,
			     valueListInt,
			     *nItems ) ;

      byte_order[ filePtr ] = Wrong_Endian ;

      if ( ierr ) LastHeaderNotFound = true;
     // printf("serial readheader_ done: %s\n", keyphrase);
    
  return ;
    }//end of it's serial
}


void
readdatablock_( int*  fileDescriptor,
                const char keyphrase[],
                void* valueArray,
                int*  nItems,
                const char  datatype[],
                const char  iotype[] ) {
  int i = *fileDescriptor;

  if(rbIONextIndex == 0)
    {
//printf("serial readdatablock %s - rank = %d\n", keyphrase, rbIOFileA[i]->rankInComm);
      int filePtr = *fileDescriptor - 1;
      FILE* fileObject;
      char junk;
    
      if ( *fileDescriptor < 1 || *fileDescriptor > (int)fileArray.size() ) {
	fprintf(stderr,"No file associated with Descriptor %d\n",*fileDescriptor);
	fprintf(stderr,"openfile_ function has to be called before \n") ;
	fprintf(stderr,"acessing the file\n ") ;
	fprintf(stderr,"fatal error: cannot continue, returning out of call\n");
	return;
      }
    
      // error check..
      // since we require that a consistant header always preceed the data block
      // let us check to see that it is actually the case.    

      if ( ! cscompare( LastHeaderKey[ filePtr ], keyphrase ) ) {
	fprintf(stderr, "Header not consistant with data block\n");
	fprintf(stderr, "Header: %s\n", LastHeaderKey[ filePtr ] );
	fprintf(stderr, "DataBlock: %s\n ", keyphrase ); 
	fprintf(stderr, "Please recheck read sequence \n");
	if( Strict_Error ) {
	  fprintf(stderr, "fatal error: cannot continue, returning out of call\n"); 
	  return;
	}
  
      }
    
      if ( LastHeaderNotFound ) return;

      fileObject = fileArray[ filePtr ];
      Wrong_Endian = byte_order[ filePtr ];

      size_t type_size = typeSize( datatype ); 
      int nUnits = *nItems;
      isBinary( iotype );
    
      if ( binary_format ) {
	fread( valueArray, type_size, nUnits, fileObject );
	fread( &junk, sizeof(char), 1 , fileObject );
	if ( Wrong_Endian ) SwapArrayByteOrder_( valueArray, type_size, nUnits );
      } else { 

	char* ts1 = StringStripper( datatype );
	if ( cscompare( "integer", ts1 ) ) {
	  for( int n=0; n < nUnits ; n++ ) 
	    fscanf(fileObject, "%d\n",(int*)((int*)valueArray+n) );  
	} else if ( cscompare( "double", ts1 ) ) {
	  for( int n=0; n < nUnits ; n++ ) 
	    fscanf(fileObject, "%lf\n",(double*)((double*)valueArray+n) );  
	}
	delete [] ts1;
      }
      //printf("serial readdatablock_ done: %s\n", keyphrase);
      return;
    }//end of it's serial
}

void
writeheaderrb_ (  int*  fileDescriptor,
		const char keyphrase[],
		void* valueArray,
		int* nItems,
		int* ndataItems,
		const char datatype[],
		const char iotype[]  ) {

  int i = *fileDescriptor;
 
  if(rbIONextIndex == 0)//it's serial lib
    {
      int filePtr = *fileDescriptor - 1;
      FILE* fileObject;


      if ( *fileDescriptor < 1 || *fileDescriptor > (int)fileArray.size() ) {
	fprintf(stderr,"No file associated with Descriptor %d\n",*fileDescriptor);
	fprintf(stderr,"openfile_ function has to be called before \n") ;
	fprintf(stderr,"acessing the file\n ") ;
	fprintf(stderr,"fatal error: cannot continue, returning out of call\n");
	return;
      }

      LastHeaderKey[ filePtr ] = const_cast< char* >( keyphrase );
      DataSize = *ndataItems;
      fileObject = fileArray[ filePtr ] ;
      size_t type_size = typeSize( datatype );
      isBinary( iotype );
      header_type[ filePtr ] = type_size;

      int _newline = ( *ndataItems > 0 ) ? sizeof( char ) : 0;
      unsigned int size_of_nextblock =
	( binary_format ) ? type_size*( *ndataItems )+ _newline : *ndataItems ;

      fprintf( fileObject, "%s : < %u > ", keyphrase, size_of_nextblock );
      for( int i = 0; i < *nItems; i++ )
	fprintf(fileObject, "%d ", *((int*)((int*)valueArray+i)));
      fprintf(fileObject, "\n");

      //printf("serial writeheader_ done: %s\n", keyphrase);////////////////////////////////////////////////////////////////
      return ;
    }
  //............................................................................................
  else //it's parallel IO
    {
      if(rbIOFileA[i]->rankInGroup != rbIOFileA[i]->writerPos) //it's a worker
	{
  
	  DataSize = *ndataItems;
	  size_t type_size = typeSize( datatype );
	  isBinary( iotype );
  
	  // int _newline = sizeof(char);
	  int _newline = 0; //don't know want to plus sizeof(char)..for"\n"?
	  unsigned int size_of_nextblock = type_size * (*ndataItems) + _newline;

	  bzero((void*)rbIOFileA[i]->datablock_header,DB_HEADER_SIZE);
  
	  sprintf( rbIOFileA[i]->datablock_header,
		   "\n%s : < %u >",
		   keyphrase,
		   size_of_nextblock );

	  for ( int j = 0; j < *nItems; j++ )
	    {
	      sprintf( rbIOFileA[i]->datablock_header,
		       "%s %d ",
		       rbIOFileA[i]->datablock_header,
		       *((int*)((int*)valueArray+j)));
	    }
	  sprintf( rbIOFileA[i]->datablock_header,
		   "%s\n ",
		   rbIOFileA[i]->datablock_header );
	  //printf("datablock_header's content is %s \n", rbIOFileA[i]->datablock_header);
          
	  //if isize is zero, phasta won't call writedatablock, so we call it here explicitly
	  if( *ndataItems == 0)
	  {
		double* placeholder;
		writedatablock_(fileDescriptor, keyphrase, (void*)placeholder, ndataItems, datatype, iotype);
 	  }	
	}
      else //it's a writer
	{
	  //do nothing here, just dummy
	  //if size is zero, phasta won't call writedatablock function, so we call it here
	  //update: this is not true!! *writers* always call both functions to make things easy
	  //if( *ndataItems == 0)
       	  //writerreceive_(fileDescriptor, datatype, iotype);
	}

    }
}

void
writedatablock_( int* fileDescriptor,
                    const char keyphrase[],
                    void* valueArray,
                    int* nItems,
                    const char datatype[],
                    const char iotype[] ) {

  int isize = *nItems;
  int  i = *fileDescriptor;

  isBinary(iotype);

  if(rbIONextIndex == 0) //it's serial lib
    {
         
      int filePtr = *fileDescriptor - 1;
    
      if ( *fileDescriptor < 1 || *fileDescriptor > (int)fileArray.size() ) {
	fprintf(stderr,"No file associated with Descriptor %d\n",*fileDescriptor);
	fprintf(stderr,"openfile_ function has to be called before \n") ;
	fprintf(stderr,"acessing the file\n ") ;
	fprintf(stderr,"fatal error: cannot continue, returning out of call\n");
	return;
      }
    
      // error check..
      // since we require that a consistant header always preceed the data block
      // let us check to see that it is actually the case.    
    
      if ( ! cscompare( LastHeaderKey[ filePtr ], keyphrase ) ) {
	fprintf(stderr, "Header not consistant with data block\n");
	fprintf(stderr, "Header: %s\n", LastHeaderKey[ filePtr ] );
	fprintf(stderr, "DataBlock: %s\n ", keyphrase ); 
	fprintf(stderr, "Please recheck write sequence \n");
	if( Strict_Error ) {
	  fprintf(stderr, "fatal error: cannot continue, returning out of call\n"); 
	  return;
	}
      }
    
      FILE* fileObject =  fileArray[ filePtr ] ;
      size_t type_size=typeSize( datatype );
      isBinary( iotype );

      if ( header_type[filePtr] != (int)type_size ) {
	fprintf(stderr,"header and datablock differ on typeof data in the block for\n"); 
	fprintf(stderr,"keyphrase : %s\n", keyphrase);
	if( Strict_Error ) {
	  fprintf(stderr,"fatal error: cannot continue, returning out of call\n" );
	  return;
	}
      }

      int nUnits = *nItems;

      if ( nUnits != DataSize ) {
	fprintf(stderr,"header and datablock differ on number of data items for\n"); 
	fprintf(stderr,"keyphrase : %s\n", keyphrase);
	if( Strict_Error ) {
	  fprintf(stderr,"fatal error: cannot continue, returning out of call\n" );
	  return;
	}
      }
 
      if ( binary_format ) {

	fwrite( valueArray, type_size, nUnits, fileObject );
	fprintf( fileObject,"\n");
        
      } else { 
        
	char* ts1 = StringStripper( datatype );
	if ( cscompare( "integer", ts1 ) ) {
	  for( int n=0; n < nUnits ; n++ ) 
	    fprintf(fileObject,"%d\n",*((int*)((int*)valueArray+n)));
	} else if ( cscompare( "double", ts1 ) ) {
	  for( int n=0; n < nUnits ; n++ ) 
	    fprintf(fileObject,"%lf\n",*((double*)((double*)valueArray+n)));
	}
	delete [] ts1;
      }
      //printf("serial writedatablock_ done: %s\n", keyphrase);//////////////////////////////////////////////////
	
      return ;
    }
  //.............................................................................
  else //it's parallel lib
    {
      MPI_Request isend_req;		
      if ( cscompare( datatype, "double") ) //if it's double
	{	
	  rbIOFileA[i]->writer_datatype = 1; //"1" stands for "double", "2" for "int"

	  if(rbIOFileA[i]->rankInGroup != rbIOFileA[i]->writerPos) //it's worker
	    {
	      //allocate another memory to copy data here and let worker go; 
  
	      //here +10 is just a magic number for extra package head, eg rank and size etc
	      rbIOFileA[i]->double_chunk[rbIOFileA[i]->field_count] = new double[DB_HEADER_SIZE/sizeof(double) + isize + 10];
	      double * field_data_send = rbIOFileA[i]->double_chunk[rbIOFileA[i]->field_count];
      
	      rbIOFileA[i]->field_count ++;
      
	      double* field_data = (double*)valueArray;

	      //#############################################
	      //augmented data, just give value just for correctness verification and debugging
	      //############################################
	      if(AUGMENT_FLAG == true)
		{
			if(isize != 0)
			{
		  field_data[0] = isize; //
		  field_data[1] =  rbIOFileA[i]->rankInComm;
		  if(DEBUG_FLAG)printf("before isend - isize for rank %d is %d\n",  rbIOFileA[i]->rankInComm, isize);//////////////////////////////////////////////////////////
		  field_data[isize -1] =  rbIOFileA[i]->rankInComm + 0.1415;  
			}
			else
			printf("!!!!ISIZE = 0, can't augment DATABLOCK part because it's empty!!!\n");
		}

	      //memcpy start timer
	      field_data_send[0] = rbIOFileA[i]->rankInGroup, field_data_send[1] = isize;//##should send rankInGroup instead of rankInComm
	      memcpy(&field_data_send[2], rbIOFileA[i]->datablock_header, DB_HEADER_SIZE); //copy db_header
	      memcpy(&field_data_send[DB_HEADER_SIZE/sizeof(double) + 2], 
		     field_data, isize * sizeof(double) );  //"+2" is counting the package header,ie rank and isize
	      //memcpy end timer

	      int position = rbIOFileA[i]->groupSize * rbIOFileA[i]-> groupRank + rbIOFileA[i]->realwriterPos;
	      MPI_Isend(field_data_send, isize + DB_HEADER_SIZE/sizeof(double) + 2, MPI_DOUBLE, position, 1, MPI_COMM_WORLD, &isend_req);	
	      if(DEBUG_FLAG)
		{
		  printf("Isend done. send size from rank %d is %d\n", rbIOFileA[i]->rankInComm, isize*sizeof(double)+DB_HEADER_SIZE);
		}
	    }
	  else if(rbIOFileA[i]->rankInGroup == rbIOFileA[i]->writerPos) 
	    //it's writer, double/int issue taken care of in writerreceive func
	    {
	      //call writerreceive function
	      //printf("###writer going into writerreceive() func\n");
	      writerreceive_(fileDescriptor, datatype, iotype);
	    }

//	 printf("parallel writedatablock_ done:%s - rank %d\n", keyphrase, rbIOFileA[i]->rankInComm);
	}
      //-----------------------------------------------------------------------------------------------------
      else if( cscompare( datatype, "integer") ) //if it's integer
	{
	  rbIOFileA[i]->writer_datatype = 2; //"1" stands for "double", "2" for "int"

	  if(rbIOFileA[i]->rankInGroup != rbIOFileA[i]->writerPos) //it's worker
	    {
	      //allocate another memory to copy data here and let worker go; 
  
	      //here +10 is just a magic number for extra package head, eg rank and size etc
	      rbIOFileA[i]->int_chunk[rbIOFileA[i]->field_count] = new int[DB_HEADER_SIZE/sizeof(int) + isize + 10];
	      int * field_data_send = rbIOFileA[i]->int_chunk[rbIOFileA[i]->field_count];
      
	      rbIOFileA[i]->field_count ++;
      
	      int* field_data = (int*)valueArray;

	      //#############################################
	      //augmented data, just give value just for correctness verification and debugging
	      //############################################
	      if(AUGMENT_FLAG == true)
		{
		  field_data[0] = isize; //
		  field_data[1] =  rbIOFileA[i]->rankInComm;
		  if(DEBUG_FLAG)printf("before isend - isize for rank %d is %d\n",  rbIOFileA[i]->rankInComm, isize);//////////////////////////////////////////////////////////
		  field_data[isize -1] =  rbIOFileA[i]->rankInComm + 23;  //
		}

	      //memcpy start timer
	      field_data_send[0] = rbIOFileA[i]->rankInGroup, field_data_send[1] = isize;//##should send rankInGroup instead of rankInComm
	      memcpy(&field_data_send[2], rbIOFileA[i]->datablock_header, DB_HEADER_SIZE); //copy db_header
	      memcpy(&field_data_send[DB_HEADER_SIZE/sizeof(int) + 2], 
		     field_data, isize * sizeof(int) );  //"+2" is counting the package header,ie rank and isize
	      //memcpy end timer

	      int position = rbIOFileA[i]->groupSize * rbIOFileA[i]-> groupRank + rbIOFileA[i]->realwriterPos;
	      MPI_Isend(field_data_send, isize + DB_HEADER_SIZE/sizeof(int) + 2, MPI_INT, position, 1, MPI_COMM_WORLD, &isend_req);	
	      if(DEBUG_FLAG)
		{
		  printf("Isend done. send size from rank %d is %d\n", rbIOFileA[i]->rankInComm, isize*sizeof(int)+DB_HEADER_SIZE);
		}
	    }
	  else if(rbIOFileA[i]->rankInGroup == rbIOFileA[i]->writerPos) 
	    //it's writer, double/int issue taken care of in writerreceive func
	    {
	      //call writerreceive function
	      writerreceive_(fileDescriptor, datatype, iotype);
	    }

	}
      else
	{
	  printf("ILLEGAL DATA TYPE!!\n");
	  //proceed with exception handling here
	}
      // ....................................................................
    }
}


void
writerreceive_(int* fileDescriptor,
	       const char datatype[],
	       const char iotype[]
	       ){
  int j = *fileDescriptor;

  MPI_Request irecv_req;
  MPI_Status irecv_status;
  int isize; //approximate isize for some message

  int *iSizeArray = new int[ rbIOFileA[j]->groupSize];      //calculate total size for big buffer write
  int *iMsgArray = new int[rbIOFileA[j]-> groupSize];  
  
  if( rbIOFileA[j]->field_count == 0)
    {
       rbIOFileA[j]->data_location = MASTER_HEADER_SIZE;
       rbIOFileA[j]->totalFileSize = MASTER_HEADER_SIZE;
      
    }
  if(DEBUG_FLAG)printf("ready to receive msg I'm rank %d, writerPos is %d \n",  rbIOFileA[j]->rankInComm,  rbIOFileA[j]->writerPos);//////////////////////////

  if( cscompare(datatype, "double") )
    {
      //##############################
      //can't help putting a timer here
      long long func_start_timer = rdtsc();

      rbIOFileA[j]->writer_datatype = 1; //"1" stands for "double", "2" for "int"

     //if it's first field, alloc big buffer for write_file func
    // if( rbIOFileA[j]->field_count == 0)
//	{
//	   rbIOFileA[j]->writer_double_buffer = new double[MAX_WRITER_BUFFER_SIZE_MB/sizeof(double)];
//	   printf("new writer_double_buffer succeed\n");
//	}
 
      //@@@@@@@@@@@@@@@@@@
      //timing stuff
      if(TIMING_FLAG)
	{
	   rbIOFileA[j]->irecv_start = new long long[ rbIOFileA[j]->groupSize];
	   rbIOFileA[j]->irecv_end = new long long[ rbIOFileA[j]->groupSize];
	   rbIOFileA[j]->irecv_time = new double[ rbIOFileA[j]->groupSize];
	   rbIOFileA[j]->irecv_start[0] = rdtsc();
	}
      //#############################


      //#################################
      //receive data packages from workers
      //note: iMsgArray[0] and iSizeArray[0] would not have value

      //################################
      int temp_first_size = ONE_MEGABYTE; //4M*8 = 32M actually
      double* temp_first = new double[temp_first_size];//use 1*8 = 8MB to receive the first msg
      
      MPI_Status recv_sta;  
      MPI_Recv(temp_first, temp_first_size, MPI_DOUBLE, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &recv_sta);
      //MPI_Wait(&irecv_req, &irecv_status);
      int iRank = (int) temp_first[0];
      iMsgArray[iRank] = 0;
      iSizeArray[iRank] = (int) temp_first[1] + DB_HEADER_SIZE/sizeof(double); //this is the total number of doubles
      isize = (int) temp_first[1];
      //check if we got the whole first package
      if(isize >= temp_first_size) 
	{
	  printf("error!!!not enough buffer for first message!\n");
	  //abort();
	}

      //###############################
      //timing
      if(TIMING_FLAG)
	{
	   rbIOFileA[j]->irecv_end[0] = rdtsc();
	   rbIOFileA[j]->irecv_time[0] = ( rbIOFileA[j]->irecv_end[0] -  rbIOFileA[j]->irecv_start[0])/700000000.0;
      
	  if(DEBUG_FLAG)
	  printf("1 T_irecv for writer %4d (with rank %4d) is %lf, from rank %4d, position %d, i = 0\n",  
		  rbIOFileA[j]->groupRank, 
		  rbIOFileA[j]->rankInComm, 
		 rbIOFileA[j]->irecv_time[0], 
		 (int) temp_first[DB_HEADER_SIZE/sizeof(double) + 3],
		 (int) temp_first[DB_HEADER_SIZE/sizeof(double) + 3] - rbIOFileA[j]->rankInComm );
	}
      //@@@@@@@@@@@@@
      //now I know approx. isize, start allocating buffers
      double *field_data_send = new double[ iSizeArray[iRank] * 3];// 2x should be enough...
      double* field_data_recv[ rbIOFileA[j]->groupSize];
      for(int i = 0 ; i < rbIOFileA[j]->groupSize; i++) field_data_recv[i] = NULL;
   
      //now moving real size first msg 
      field_data_recv[iRank] = new double[iSizeArray[iRank]];
    
      if(DEBUG_FLAG)
	printf("first msg!!allocating field_data_recv[%d] in group %d addr %x, isize is %d\n", iRank,rbIOFileA[j]->groupRank, field_data_recv[iRank], iSizeArray[iRank]);

      memcpy(&field_data_recv[iRank][0], &temp_first[2], iSizeArray[iRank] * sizeof(double)); //rememer iSizeArray has total no.
      delete [] temp_first;//no need to keep it, we got approx. isize now
      
      //now processs rest messages
      for( int i = 1; i <  rbIOFileA[j]->groupSize -1; i ++)
	{
	  //################################
	  if(TIMING_FLAG)
	    rbIOFileA[j]->irecv_start[i] = rdtsc();
	  //#################################

	  MPI_Recv(field_data_send, 
		    (DB_HEADER_SIZE/sizeof(double) + 2) + isize *3, 
		    MPI_DOUBLE, 
		    MPI_ANY_SOURCE, 
		    1, 
		    MPI_COMM_WORLD, 
		    &recv_sta);
	  //MPI_Wait(&irecv_req, &irecv_status);
	  iRank = (int) field_data_send[0];
	  iMsgArray[iRank] = i;
	  iSizeArray[iRank] = (int)  field_data_send[1] + DB_HEADER_SIZE/sizeof(double);

	  field_data_recv[iRank] = new double[iSizeArray[iRank]];

	  if(DEBUG_FLAG)
	    printf("allocating field_data_recv[%d] in group %d addr %x, isize is %d\n", 
		   iRank,   rbIOFileA[j]->groupRank, field_data_recv[iRank], iSizeArray[iRank]);
	  
	  //####################################
	  if(TIMING_FLAG)
	    {
	      rbIOFileA[j]->irecv_end[i] = rdtsc();
	      rbIOFileA[j]->irecv_time[i] = (rbIOFileA[j]->irecv_end[i] - rbIOFileA[j]->irecv_start[i])/700000000.0;
      
	  if(DEBUG_FLAG)
	  printf("1 T_irecv for writer %4d (with rank %4d) is %lf, from rank %4d, position %d, i = %d\n", 
		 rbIOFileA[j]->groupRank,  
		 rbIOFileA[j]->rankInComm, 
		 rbIOFileA[j]->irecv_time[i], 
		 (int) field_data_send[DB_HEADER_SIZE/sizeof(double) + 3],
		 (int) field_data_send[DB_HEADER_SIZE/sizeof(double) + 3] -  rbIOFileA[j]->rankInComm , 
		 i);
	    }
	  //######################################

	  memcpy(&field_data_recv[iRank][0], &field_data_send[2], iSizeArray[iRank] * sizeof(double) );
	
	  if(DEBUG_FLAG)
	    printf("receive size is %d for rank %d in group %d\n", 
		   iSizeArray[iRank]*sizeof(double), iRank, rbIOFileA[j]->groupRank);
	}//end of for loop for irecv

	if(DEBUG_FLAG)
	printf("done with the loop\n");

	/**************************************/   
      //###############################################
      //compute total time of irecvs
      double total_irecv_time = 0;
      for( int i = 0; i < rbIOFileA[j]->groupSize -1; i++)
	total_irecv_time += rbIOFileA[j]->irecv_time[i];
      //printf("total T_irecv for writer %4d (with rank %4d) is %lf\n", groupRank, rankInComm, total_irecv_time);

      //###############################
      //fill information for master header
      //###############################
      
      //update before MPI_write and flush before MPI_File_close
		
      //################################
      //fill the offset table to offset_table[]
      //###############################
   
      long long* offset_table = new long long[rbIOFileA[j]->groupSize];
      long long numItems = 0;
      for( int i = 1; i < rbIOFileA[j]->groupSize ; i ++) //same as irecv
	{
	  numItems += iSizeArray[i] ;
	  //printf("iSizeArray[%d] is %d for rank \n", i, iSizeArray[i]*sizeof(double));
	  offset_table[i] =  ( rbIOFileA[j]->data_location + (numItems - iSizeArray[i])*sizeof(double));
	  //if(groupRank ==0)printf("offset table[%d] is %lld \n", i, offset_table[i]);
	}
      //cout << "total data is "<< numItems * sizeof(double)/1000000<< "M bytes" <<endl;
      //printf("numItems is %lld total data is %lld M bytes \n", numItems, numItems * sizeof(double) / 1000000);
        	
      
      rbIOFileA[j]->thisFieldSize = numItems * sizeof(double);//in bytes
      int allDataSize = rbIOFileA[j]->thisFieldSize / sizeof(double) + MASTER_HEADER_SIZE/sizeof(double) ;
      

	if(DEBUG_FLAG)
        printf("allData new to be done. , Size is %d \n",  allDataSize);

	double* allData = new double[allDataSize ];	
      if(DEBUG_FLAG)
	printf("allData new done. allData addr is %x, Size is %d \n", allData, allDataSize);

      //#################################
      //assemble data into a big buffer
      //#################################
      // after version info comes offset table
      int temp_offset_location = VERSION_INFO_HEADER_SIZE/sizeof(double);
      memcpy(&allData[temp_offset_location], 
	     &offset_table[1], 
	     ( rbIOFileA[j]->groupSize-1)*sizeof(long long));//from offset_table[1] not [0] because 0 is writerPos

	
      //finally the real data
      long long Index = 0;
      for( int i = 1; i < rbIOFileA[j]->groupSize ; i++)
	{
	  //for rank i, msg# is i, size info is in iSizeArray[i]
	  memcpy(&(allData[((MASTER_HEADER_SIZE/sizeof(double)) + Index)]), field_data_recv[ i ], iSizeArray[i] * sizeof(double) );
	  //memcpy is measured in bytes not number of doubles
	  Index += iSizeArray[i];//Index is measured in number of doubles, i.e. it's allData[Index].
	  delete [] field_data_recv[i];
	  if(DEBUG_FLAG) printf("memcpy from field_data_recv[%d] to allData done\n", i);
	}
      int ilast = rbIOFileA[j]->groupSize;
	
      //##################################
      //update master header info and data_location info etc
      //##################################
      
      //update field_table in version info
      //decide to parse out field_name from messge thus save user the trouble of giving field_name in arguments
      
      char parse_field_name[MAX_FIELDS_NAME_LENGTH];
      memcpy(parse_field_name, &field_data_send[2], MAX_FIELDS_NAME_LENGTH);
	if(DEBUG_FLAG) {for (int i = 0; i < MAX_FIELDS_NAME_LENGTH; i++)printf("%c ",parse_field_name[i]);}
      char* parsed = strtok(parse_field_name, "@");
      if(DEBUG_FLAG) printf("here parsed is %s\n", parsed+1);
      //use "parsed+1" coz' to get rid of the "\n" at the beginning of every db header
      //refer to phastaio spec. fail to +1 would cause reading failure
      sprintf(  rbIOFileA[j]->field_table[ rbIOFileA[j]->field_count],"%s:", parsed+1);
      if(DEBUG_FLAG) printf("field parsing done - %s\n", rbIOFileA[j]->field_table[ rbIOFileA[j]->field_count ]);

      //update offset_table
      memcpy( rbIOFileA[j]->my_offset_table[ rbIOFileA[j]->field_count], &offset_table[1], ( rbIOFileA[j]->groupSize-1)*sizeof(long long));
      
      if(DEBUG_FLAG)
	{
	  for( int i = 0; i <  rbIOFileA[j]->groupSize -1; i++) printf(" %lld == %lld \n",  rbIOFileA[j]->my_offset_table[ rbIOFileA[j]->field_count][i], offset_table[1+i]);
	}

	if( rbIOFileA[j]->field_count == 0)
        {
           rbIOFileA[j]->writer_double_buffer = new double[MAX_WRITER_BUFFER_SIZE_MB/sizeof(double)];
           //printf("new writer_double_buffer succeed\n");
        }

      //update field number we got
       rbIOFileA[j]->field_count ++;

      //update file size
       rbIOFileA[j]->totalFileSize +=  rbIOFileA[j]->thisFieldSize;

      //copy to global array
      if(DEBUG_FLAG)
	{
	  printf("data_location is %lld, thisFieldSize is %lld\n",  rbIOFileA[j]->data_location,  rbIOFileA[j]->thisFieldSize);
	}
      memcpy(&( rbIOFileA[j]->writer_double_buffer[( (rbIOFileA[j]->data_location ) / sizeof(double) )]), &(allData[ (MASTER_HEADER_SIZE/sizeof(double)) ]),  rbIOFileA[j]->thisFieldSize);
     
      //printf("write total file size of %lld\n", totalFileSize);

      rbIOFileA[j]->data_location =  rbIOFileA[j]->totalFileSize;
      
      //#################################
      //clean things up, free memory
      //#################################
         
       for( int i = 1 ; i < rbIOFileA[j]->groupSize ; i++)
	{
	  //free(field_data_recv[i]);
	  if(field_data_recv[i] != NULL)
	    {
	     // delete [] field_data_recv[i];
	      field_data_recv[i] = NULL;
	    }
	}
       
       delete [] iMsgArray;
       delete [] offset_table;


      delete []  allData;
      //allData = NULL;
      delete []  iSizeArray;
      //iSizeArray = NULL;
      delete [] field_data_send;
      //field_data_send = NULL;
       //#################################
      //compute stage timings
      //#################################
      long long func_end_timer = rdtsc();
      double func_time = (func_end_timer - func_start_timer)/700000000.0;
      
      if(DEBUG_FLAG)
      printf("total writerreceive func time for writer %4d (with rank %4d) is %lf, misc time is %lf\n",  rbIOFileA[j]->groupRank,  rbIOFileA[j]->rankInComm, func_time, func_time - total_irecv_time);
   
      return;	
    }//end of if datatype is double 
  //----------------------------------------------------------------------------------
  else if ( cscompare( datatype, "integer") )
    {

      //##############################
      //can't help putting a timer here
      long long func_start_timer = rdtsc();
      
      rbIOFileA[j]->writer_datatype = 2; //"1" stands for "double", "2" for "int"

      //if it's first field, alloc big buffer for write_file func
      if( rbIOFileA[j]->field_count == 0)
	{
	   rbIOFileA[j]->writer_int_buffer = new int[MAX_WRITER_BUFFER_SIZE_MB/sizeof(int)];
	}

       //@@@@@@@@@@@@@@@@@@
      //timing stuff
      if(TIMING_FLAG)
	{
	   rbIOFileA[j]->irecv_start = new long long[ rbIOFileA[j]->groupSize];
	   rbIOFileA[j]->irecv_end = new long long[ rbIOFileA[j]->groupSize];
	   rbIOFileA[j]->irecv_time = new double[ rbIOFileA[j]->groupSize];
	   rbIOFileA[j]->irecv_start[0] = rdtsc();
	}
      //#############################

      //#################################
      //receive data packages from workers
      //note: iMsgArray[0] and iSizeArray[0] would not have value

      //################################
      int temp_first_size = ONE_MEGABYTE; //4M*8 = 32M actually
      int* temp_first = new int[temp_first_size];//use 1*8 = 8MB to receive the first msg
            
      MPI_Irecv(temp_first, temp_first_size, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &irecv_req);
      MPI_Wait(&irecv_req, &irecv_status);
      int iRank = (int) temp_first[0];
      iMsgArray[iRank] = 0;
      iSizeArray[iRank] = (int) temp_first[1] + DB_HEADER_SIZE/sizeof(int); //this is the total number of doubles
      isize = (int) temp_first[1];
      //check if we got the whole first package
      if(isize >= temp_first_size) 
	{
	  printf("error!!!not enough buffer for first message!\n");
	  //abort();
	}

      //###############################
      //timing
      if(TIMING_FLAG)
	{
	   rbIOFileA[j]->irecv_end[0] = rdtsc();
	   rbIOFileA[j]->irecv_time[0] = ( rbIOFileA[j]->irecv_end[0] -  rbIOFileA[j]->irecv_start[0])/700000000.0;
      
	  if(DEBUG_FLAG)
	  printf("1 T_irecv for writer %4d (with rank %4d) is %lf, from rank %4d, position %d, i = 0\n",  
		  rbIOFileA[j]->groupRank, 
		  rbIOFileA[j]->rankInComm, 
		 rbIOFileA[j]->irecv_time[0], 
		 (int) temp_first[DB_HEADER_SIZE/sizeof(int) + 3],
		 (int) temp_first[DB_HEADER_SIZE/sizeof(int) + 3] - rbIOFileA[j]->rankInComm );
	}
      //@@@@@@@@@@@@@

      //now I know approx. isize, start allocating buffers
      int *field_data_send = new int[ isize * 3];// 2x should be enough...
    
      int * field_data_recv[ rbIOFileA[j]->groupSize];
      for(int i = 0 ; i < rbIOFileA[j]->groupSize; i++) field_data_recv[i] = NULL;

      //now moving real size first msg 
      field_data_recv[iRank] = new int[iSizeArray[iRank]];
    
      if(DEBUG_FLAG)
	printf("first msg!!allocating field_data_recv[%d] in group %d addr %x, isize is %d\n", iRank,rbIOFileA[j]->groupRank, field_data_recv[iRank], iSizeArray[iRank]);

      memcpy(&field_data_recv[iRank][0], &temp_first[2], iSizeArray[iRank] * sizeof(int)); //rememer iSizeArray has total no.
      delete [] temp_first;//no need to keep it, we got approx. isize now

    //now processs rest messages
      for( int i = 1; i <  rbIOFileA[j]->groupSize -1; i ++)
	{
	  //################################
	  if(TIMING_FLAG)
	    rbIOFileA[j]->irecv_start[i] = rdtsc();
	  //#################################

	  MPI_Irecv(field_data_send, 
		    isize *3, 
		    MPI_INT, 
		    MPI_ANY_SOURCE, 
		    1, 
		    MPI_COMM_WORLD, 
		    &irecv_req);
	  MPI_Wait(&irecv_req, &irecv_status);
	  iRank = (int) field_data_send[0];
	  iMsgArray[iRank] = i;
	  iSizeArray[iRank] = (int)  field_data_send[1] + DB_HEADER_SIZE/sizeof(int);

	  field_data_recv[iRank] = new int[iSizeArray[iRank]];

	  if(DEBUG_FLAG)
	    printf("allocating field_data_recv[%d] in group %d addr %x, isize is %d\n", 
		   iRank, rbIOFileA[j]->groupRank, field_data_recv[iRank], iSizeArray[iRank]);
	  
	  //####################################
	  if(TIMING_FLAG)
	    {
	      rbIOFileA[j]->irecv_end[i] = rdtsc();
	      rbIOFileA[j]->irecv_time[i] = (rbIOFileA[j]->irecv_end[i] - rbIOFileA[j]->irecv_start[i])/700000000.0;
      
	  if(DEBUG_FLAG)
	    printf("1 T_irecv for writer %4d (with rank %4d) is %lf, from rank %4d, position %d, i = %d\n", 
		 rbIOFileA[j]->groupRank,  
		 rbIOFileA[j]->rankInComm, 
		 rbIOFileA[j]->irecv_time[i], 
		 (int) field_data_send[DB_HEADER_SIZE/sizeof(int) + 3],
		 (int) field_data_send[DB_HEADER_SIZE/sizeof(int) + 3] -  rbIOFileA[j]->rankInComm , 
		 i);
	    }
	  //######################################

	  memcpy(&field_data_recv[iRank][0], &field_data_send[2], iSizeArray[iRank] * sizeof(int) );
	
	  if(DEBUG_FLAG)
	    printf("receive size is %d for rank %d in group %d\n", 
		   iSizeArray[iRank]*sizeof(int), iRank, rbIOFileA[j]->groupRank);
	}//end of for loop for irecv

  
      //###############################################
      //compute total time of irecvs
      double total_irecv_time = 0;
      for( int i = 0; i < rbIOFileA[j]->groupSize -1; i++)
	total_irecv_time += rbIOFileA[j]->irecv_time[i];
      //printf("total T_irecv for writer %4d (with rank %4d) is %lf\n", groupRank, rankInComm, total_irecv_time);

      //###############################
      //fill information for master header
      //###############################
      
      //update before MPI_write and flush before MPI_File_close
		
      //################################
      //fill the offset table to offset_table[]
      //###############################

      long long* offset_table = new long long[rbIOFileA[j]->groupSize];
      long long numItems = 0;
      for( int i = 1; i < rbIOFileA[j]->groupSize ; i ++) //same as irecv
	{
	  numItems += iSizeArray[i] ;
	  //printf("iSizeArray[%d] is %d for rank \n", i, iSizeArray[i]*sizeof(double));
	  offset_table[i] =  ( rbIOFileA[j]->data_location + (numItems - iSizeArray[i])*sizeof(int));
	  //if(groupRank ==0)printf("offset table[%d] is %lld \n", i, offset_table[i]);
	}
      //cout << "total data is "<< numItems * sizeof(int)/1000000<< "M bytes" <<endl;
      //printf("numItems is %lld total data is %lld M bytes \n", numItems, numItems * sizeof(int) / 1000000);
        
     
      rbIOFileA[j]->thisFieldSize = numItems * sizeof(int);//in bytes
      int allDataSize = rbIOFileA[j]->thisFieldSize / sizeof(int) + MASTER_HEADER_SIZE/sizeof(int) ;
      int* allData = new int[allDataSize ];	
      //printf("allData addr is %x \n", allData);

      //#################################
      //assemble data into a big buffer
      //#################################
      // after version info comes offset table
      int temp_offset_location = VERSION_INFO_HEADER_SIZE/sizeof(int);
      memcpy(&allData[temp_offset_location], 
	     &offset_table[1], 
	     ( rbIOFileA[j]->groupSize-1)*sizeof(long long));//from offset_table[1] not [0] because 0 is writerPos

      //finally the real data
      long long Index = 0;
      for( int i = 1; i < rbIOFileA[j]->groupSize ; i++)
	{
	  //for rank i, msg# is i, size info is in iSizeArray[i]
	  memcpy(&(allData[((MASTER_HEADER_SIZE/sizeof(int)) + Index)]), field_data_recv[ i ], iSizeArray[i] * sizeof(int) );
	  //memcpy is measured in bytes not number of doubles
	  Index += iSizeArray[i];//Index is measured in number of doubles, i.e. it's allData[Index].
	}
      int ilast = rbIOFileA[j]->groupSize;

      //##################################
      //update master header info and data_location info etc
      //##################################
      
      //update field_table in version info
      //decide to parse out field_name from messge thus save user the trouble of giving field_name in arguments
      
      char parse_field_name[MAX_FIELDS_NAME_LENGTH];
      memcpy(parse_field_name, &field_data_send[2], MAX_FIELDS_NAME_LENGTH);
      char* parsed = strtok(parse_field_name, "@");

      //use "parsed+1" coz' to get rid of the "\n" at the beginning of every db header
      //refer to phastaio spec. fail to +1 would cause reading failure
      sprintf(  rbIOFileA[j]->field_table[ rbIOFileA[j]->field_count],"%s:", parsed+1);

      //update offset_table
      memcpy( rbIOFileA[j]->my_offset_table[ rbIOFileA[j]->field_count], &offset_table[1], ( rbIOFileA[j]->groupSize-1)*sizeof(long long));

     
      if(DEBUG_FLAG)
	{
	  for( int i = 0; i <  rbIOFileA[j]->groupSize -1; i++) printf(" %lld == %lld \n",  rbIOFileA[j]->my_offset_table[ rbIOFileA[j]->field_count][i], offset_table[1+i]);
	}

      //update field number we got
       rbIOFileA[j]->field_count ++;

      //update file size
       rbIOFileA[j]->totalFileSize +=  rbIOFileA[j]->thisFieldSize;

      //copy to global array
      if(DEBUG_FLAG)
	{
	  printf("data_location is %lld, thisFieldSize is %lld\n",  rbIOFileA[j]->data_location,  rbIOFileA[j]->thisFieldSize);
	}
      memcpy(&( rbIOFileA[j]->writer_double_buffer[( (rbIOFileA[j]->data_location ) / sizeof(int) )]), &(allData[ (MASTER_HEADER_SIZE/sizeof(int)) ]),  rbIOFileA[j]->thisFieldSize);


     //printf("write total file size of %lld\n", totalFileSize);

      rbIOFileA[j]->data_location =  rbIOFileA[j]->totalFileSize;
      
      //#################################
      //clean things up, free memory
      //#################################

       for( int i = 1 ; i < rbIOFileA[j]->groupSize ; i++)
	{
	  //free(field_data_recv[i]);
	  if(field_data_recv[i] != NULL)
	    {
	      delete [] field_data_recv[i];
	      field_data_recv[i] = NULL;
	    }
	}
      
       delete [] iMsgArray;
       delete [] offset_table;

      delete []  allData;
      //allData = NULL;
      delete []  iSizeArray;
      //iSizeArray = NULL;
      delete [] field_data_send;
      //field_data_send = NULL;
       //#################################
      //compute stage timings
      //#################################
      long long func_end_timer = rdtsc();
      double func_time = (func_end_timer - func_start_timer)/700000000.0;
      
      if(DEBUG_FLAG)
      printf("total writerreceive func time for writer %4d (with rank %4d) is %lf, misc time is %lf\n",  
	     rbIOFileA[j]->groupRank,  rbIOFileA[j]->rankInComm, func_time, func_time - total_irecv_time);
     
      return;	
    }//end of if datatype is int 
  else
    {
      printf("wrong datatype!! - has to be either double or integer!!\n");
    }
}



void 
writestring_( int* fileDescriptor,
              const char string[] ) {
    
  int filePtr = *fileDescriptor - 1;
  FILE* fileObject = fileArray[filePtr] ;
  fprintf(fileObject,"%s",string );
  return;
}

void
Gather_Headers( int* fileDescriptor,
                vector< string >& headers ) {

  FILE* fileObject;
  char Line[1024];

  fileObject = fileArray[ (*fileDescriptor)-1 ];

  while( !feof(fileObject) ) {
    fgets( Line, 1024, fileObject);
    if ( Line[0] == '#' ) {
      headers.push_back( Line );
    } else { 
      break; 
    }
  }
  rewind( fileObject );
  clearerr( fileObject );
}
void
isWrong( void ) { (Wrong_Endian) ? fprintf(stdout,"YES\n"): fprintf(stdout,"NO\n") ; }

void 
togglestrictmode_( void ) { Strict_Error = !Strict_Error; }

int
isLittleEndian_( void ) { 
  // this function returns a 1 if the current running architecture is
  // LittleEndian Byte Ordered, else it returns a zero 

  union  {
    long a;
    char c[sizeof( long )];
  } endianUnion;

  endianUnion.a = 1 ;

  if ( endianUnion.c[sizeof(long)-1] != 1 ) return 1 ;
  else return 0;
}

