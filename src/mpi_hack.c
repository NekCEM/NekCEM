#include <unistd.h>


void gnu_getcwd (char *buffer, int *fsize)
{
  size_t size = *fsize;

  getcwd (buffer, size);
  printf("pwd=%s\n",buffer);
  flush_io();
}


void  mpi_init(void *ierr)
{
  MPI_INIT(ierr);
}

void mpi_irecv(void *work, void *len, void *HPI_BYTE, void *rgtnbr, void *jtype,
	       void *HPI_COMM_WORLD, void *msg, void *ierr)
{
  MPI_IRECV(work,len,HPI_BYTE,rgtnbr,jtype,HPI_COMM_WORLD,msg,ierr);
}

void mpi_send (void *wrk ,void *len,void *HPI_BYTE,void *lftnbr,void *itype,
	       void *HPI_COMM_WORLD,void *ierr)
{
  MPI_SEND (wrk ,len,HPI_BYTE,lftnbr,itype ,HPI_COMM_WORLD,ierr);
}

void mpi_wait(void *msg,void *istat,void *ierr)
{
  MPI_WAIT(msg,istat,ierr);
}

void mpi_allreduce(void *x,void *work,void *n,void *HPI_DOUBLE_PRECISION,void *hpi_sum,
		   void *hpi_comm_world,void *ierr)
{
  MPI_ALLREDUCE(x,work,n,HPI_DOUBLE_PRECISION,hpi_sum,hpi_comm_world,ierr);
}

void mpi_iprobe(void *hpi_any_source,void *mtype,void *HPI_COMM_WORLD,void *iflag,
		void *status,void *ierr)
{
  MPI_IPROBE(hpi_any_source,mtype,HPI_COMM_WORLD,iflag,status,ierr);
}


void mpi_recv(void *buf,void *len,void *HPI_BYTE,void *jnid,void *mtype,
	      void *HPI_COMM_WORLD,void *status,void *ierr)
{
  MPI_RECV(buf,len,HPI_BYTE,jnid,mtype,HPI_COMM_WORLD,status,ierr);
}

void mpi_comm_size(void *HPI_COMM_WORLD,void *numprocs ,void * ierr)
{
  MPI_COMM_SIZE(HPI_COMM_WORLD, numprocs , ierr);
}


void mpi_comm_rank(void *HPI_COMM_WORLD,void *myid,void *ierr)
{
  MPI_COMM_RANK(HPI_COMM_WORLD, myid, ierr);
}

/*
fortran mpi_wtime verions is double precision
C       mpi_wtime verions is double
entire code is real*4 and float (i.e. 4 byte reals)
*/
/*
void hpi_wtime(float *rval)
{
  double MPI_Wtime();

  *rval = (float) MPI_Wtime();
}
*/

double mpi_wtime()
{
  double MPI_Wtime();

  return(MPI_Wtime());
}


void mpi_bcast(void *buf,void *len,void *HPI_BYTE,void *root,void *HPI_COMM_WORLD,void *ierr)
{
  MPI_BCAST(buf,len,HPI_BYTE,root,HPI_COMM_WORLD,ierr);
}

void mpi_finalize(void)
{
  MPI_FINALIZE();
}


void mpi_comm_group(void *HPI_COMM_WORLD, void *myid, void *ierr)
{
  MPI_COMM_GROUP(HPI_COMM_WORLD,myid,ierr);
}


void mpi_comm_create(void *HPI_COMM_WORLD, void *nekgroup, void *nekcomm,  void *ierr)
{
  MPI_COMM_CREATE(HPI_COMM_WORLD,nekgroup,nekcomm,ierr);
}


void mpi_group_free(void *nekgroup, void *ierr)
{
  MPI_GROUP_FREE(nekgroup,ierr);
}

void mpi_isend(void *x, void *len, void *mpi_byte, void *jnid, void *msgtag, void *nekcomm, void *imsg, void *ierr)
{
  MPI_ISEND(x,len,mpi_byte,jnid,msgtag,nekcomm,imsg,ierr);
}
