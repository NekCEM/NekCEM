/* compile-time settings:

   FORTRAN naming convention
     default      cpgs_setup, etc.
     -DUPCASE     CPGS_SETUP, etc.
     -DUNDERSCORE cpgs_setup_, etc.

   -DMPI        parallel version (sequential otherwise)
   -DDYNAMIC    exchange twice as many messages in crystal router to
                avoid crashing b/c of insufficient buffer size
   
   -DINITIAL_BUFFER_SIZE=expression
      arithmetic expression controlling the initial buffer size for the crystal
      router; when DYNAMIC is not defined, this needs to be large enough to
      hold the intermediate messages during all stages of the crystal router
      
      variables that can be used in expression include
         num   - the number of processors
         vn    - the length of the global index array
         count - the number of non-zero entries in the global index array

      note that the buffer will be expanded to hold out-going messages
      for the crystal router regardless, and if DYNAMIC is defined,
      the crystal router will be able to expand the buffer to accomodate
      intermediate stages and the incoming messages so that
      -DDYNAMIC -DINITIAL_BUFFER_SIZE=0 would be perfectly reasonable
*/

/* defaults for INITIAL_BUFFER_SIZE */
#ifndef INITIAL_BUFFER_SIZE
#  ifndef DYNAMIC
#    define INITIAL_BUFFER_SIZE 2*(3*num+vn*9)
#  else
#    define INITIAL_BUFFER_SIZE 0
#  endif
#endif

/* FORTRAN usage:

   call cpgs_setup(global_index_array, n, max_vec_size,
                   mpi_communicator, nproc, handle)
     integer n, max_vec_size, handle
     integer global_index_array(1:n)
           ? mpi_communicator
     integer nproc
     
     it is a run-time error if MPI_Comm_size gives a value other than nproc,
     or if nproc != 1 and not compiled with MPI support (-DMPI)
     (this second check is the reason for passing the redundant nproc)

   call cpgs_op(handle, u, op)
     integer handle, op : 1-add, 2-multiply, 3-min, 4-max
     real    u(1:n) - same layout as global_index_array provided to cpgs_setup

   call cpgs_op_vec(handle, u, d, op)
     integer handle, op : 1-add, 2-multiply, 3-min, 4-max
     integer d
     real    u(1:d, 1:n) - vector components for each node stored together

   call cpgs_op_many(handle, u, v, w, d, op)
     integer handle, op : 1-add, 2-multiply, 3-min, 4-max
     integer d = 1, 2, or 3
     real    u(1:n), v(1:n), w(1:n)
     
     same effect as: call cpgs_op(handle, u, op)
                     if(d.gt.1) call cpgs_op(handle, v, op)
                     if(d.gt.2) call cpgs_op(handle, w, op)
     with possibly some savings as fewer messages are exchanged
   
   call cpgs_free(handle)
     integer handle
*/

#include <stdio.h>   /* for printf */
#include <stdarg.h>  /* to be able to implement fail */
#include <stdlib.h>  /* for malloc, calloc, realloc, free */
#include <string.h>  /* for memcpy */
#ifdef MPI
#  include <mpi.h>
#endif

#ifdef UPCASE
#  define FORTRAN_NAME(low,up) up
#else
#ifdef UNDERSCORE
#  define FORTRAN_NAME(low,up) low##_
#else
#  define FORTRAN_NAME(low,up) low
#endif
#endif

#define FORTRAN_SETUP_NAME   FORTRAN_NAME(cpgs_setup  ,CPGS_SETUP  )
#define FORTRAN_OP_NAME      FORTRAN_NAME(cpgs_op     ,CPGS_OP     )
#define FORTRAN_OP_VEC_NAME  FORTRAN_NAME(cpgs_op_vec ,CPGS_OP_VEC )
#define FORTRAN_OP_MANY_NAME FORTRAN_NAME(cpgs_op_many,CPGS_OP_MANY)
#define FORTRAN_FREE_NAME    FORTRAN_NAME(cpgs_free   ,CPGS_FREE   )

#define OP_ADD 1
#define OP_MUL 2
#define OP_MIN 3
#define OP_MAX 4

/*--------------------------------------------------------------------------
   Error Reporting
   Memory Allocation Wrappers to Catch Out-of-memory
  --------------------------------------------------------------------------*/

static void fail(const char *fmt, ...)
{
  va_list ap;
  va_start(ap, fmt);
  vfprintf(stderr, fmt, ap);
  va_end(ap);
  exit(1);
}

static void failwith(const char *string)
{
  fail("%s\n",string);
}

static void *smalloc(size_t size)
{
  void *res = malloc(size);
  if(!res) fail("gs2.c: allocation of %d bytes failed\n",(int)size);
  return res;
}

static void *scalloc(size_t nmemb, size_t size)
{
  void *res = calloc(nmemb, size);
  if(!res) fail("gs2.c: allocation of %d bytes failed\n",(int)size*nmemb);
  return res;
}

static void *srealloc(void *ptr, size_t size)
{
  void *res = realloc(ptr, size);
  if(size && !res) fail("gs2.c: allocation of %d bytes failed\n",(int)size);
  return res;
}

/*--------------------------------------------------------------------------
   Local Execution Phases
  --------------------------------------------------------------------------*/

static void local_condense(double *u, int op, const int *cm)
{
  int i,j;
  switch(op) {
    case OP_ADD: while((i=*cm++) != -1) while((j=*cm++) != -1)
                   u[i]+=u[j];
      break;
    case OP_MUL: while((i=*cm++) != -1) while((j=*cm++) != -1)
                   u[i]*=u[j];
      break;
    case OP_MIN: while((i=*cm++) != -1) while((j=*cm++) != -1)
                   if(u[j]<u[i]) u[i]=u[j];
      break;
    case OP_MAX: while((i=*cm++) != -1) while((j=*cm++) != -1)
                   if(u[j]>u[i]) u[i]=u[j];
      break;
  }
}

static void local_uncondense(double *u, const int *cm)
{
  int i,j;
  while((i=*cm++) != -1) while((j=*cm++) != -1) u[j]=u[i];
}

static void local_condense_vec(double *u, int n, int op, const int *cm)
{
  int i,j,k;
  double *ui, *uii, *uj;
  switch(op) {
    case OP_ADD: while((i=*cm++) != -1) { uii=u+n*i; while((j=*cm++) != -1) {
                   ui=uii, uj=u+n*j;
                   for(k=n;k;--k) *ui++ += *uj++;
      } } break;
    case OP_MUL: while((i=*cm++) != -1) { uii=u+n*i; while((j=*cm++) != -1) {
                   ui=uii, uj=u+n*j;
                   for(k=n;k;--k) *ui++ *= *uj++;
      } } break;
    case OP_MIN: while((i=*cm++) != -1) { uii=u+n*i; while((j=*cm++) != -1) {
                   ui=uii, uj=u+n*j;
                   for(k=n;k;--k) { if (*uj < *ui) *ui = *uj; ++ui, ++uj; }
      } } break;
    case OP_MAX: while((i=*cm++) != -1) { uii=u+n*i; while((j=*cm++) != -1) {
                   ui=uii, uj=u+n*j;
                   for(k=n;k;--k) { if (*uj > *ui) *ui = *uj; ++ui, ++uj; }
      } } break;
  }
}

static void local_uncondense_vec(double *u, int n, const int *cm)
{
  int i,j;
  size_t size=n*sizeof(double);
  while((i=*cm++) != -1) while((j=*cm++) != -1) memcpy(u+n*j,u+n*i,size);
}

/*--------------------------------------------------------------------------
   Non-local Execution Phases
  --------------------------------------------------------------------------*/

#ifdef MPI
typedef struct {
  int np;            /* number of processors to communicate with          */
  int *target;       /* int target[np]: array of processor ids to comm w/ */
  int *nshared;      /* nshared[i] = number of points shared w/ target[i] */
  int *sh_ind;       /* list of shared point indices                      */
  MPI_Request *reqs; /* pre-allocated for MPI calls                       */
  double *buf;       /* pre-allocated buffer to receive data              */
  int buflen;        /* size of buf                                       */
  int maxv;          /* maximum vector size                               */
} nonlocal_transfer_info;

static void nonlocal_transfer(double *u, int op,
                              const nonlocal_transfer_info *info, MPI_Comm comm)
{
  MPI_Status status;
  int np = info->np, i;
  MPI_Request *reqs = info->reqs;
  int *targ = info->target;
  int *nshared = info->nshared;
  int *sh_ind = info->sh_ind;
  int id;
  double *buf = info->buf, *start;
  MPI_Comm_rank(comm,&id);
  for(i=0;i<np;++i) {
    int c = nshared[i];
    start = buf;
    for(;c;--c) *buf++ = u[*sh_ind++];
    MPI_Isend(start,nshared[i],MPI_DOUBLE,targ[i],id,comm,reqs++);
  }
  start = buf;
  for(i=0;i<np;++i) {
    MPI_Irecv(start,nshared[i],MPI_DOUBLE,targ[i],targ[i],comm,reqs++);
    start+=nshared[i];
  }
  for(reqs=info->reqs,i=np*2;i;--i) MPI_Wait(reqs++,&status);
  sh_ind = info->sh_ind;
  switch(op) {
    case OP_ADD: for(i=0;i<np;++i) { int c; for(c=nshared[i];c;--c)
        u[*sh_ind++] += *buf++;
      } break;
    case OP_MUL: for(i=0;i<np;++i) { int c; for(c=nshared[i];c;--c)
        u[*sh_ind++] *= *buf++;
      } break;
    case OP_MIN: for(i=0;i<np;++i) { int c; for(c=nshared[i];c;--c) {
        if(*buf < u[*sh_ind]) u[*sh_ind] = *buf;
        ++sh_ind, ++buf;
      } } break;
    case OP_MAX: for(i=0;i<np;++i) { int c; for(c=nshared[i];c;--c) {
        if(*buf > u[*sh_ind]) u[*sh_ind] = *buf;
        ++sh_ind, ++buf;
      } } break;
  }
}

static void nonlocal_transfer_vec(
              double *u, int n, int op,
              const nonlocal_transfer_info *info, MPI_Comm comm)
{
  MPI_Status status;
  int np = info->np, i;
  MPI_Request *reqs = info->reqs;
  int *targ = info->target;
  int *nshared = info->nshared;
  int *sh_ind = info->sh_ind;
  int id;
  double *buf = info->buf, *start;
  int size = n*sizeof(double);
  MPI_Comm_rank(comm,&id);
  for(i=0;i<np;++i) {
    int ns=nshared[i]; int c=ns;
    start = buf;
    for(;c;--c) memcpy(buf,u+n*(*sh_ind++),size), buf+=n;
    MPI_Isend(start,ns*n,MPI_DOUBLE,targ[i],id,comm,reqs++);
  }
  start = buf;
  for(i=0;i<np;++i) {
    int nsn=n*nshared[i];
    MPI_Irecv(start,nsn,MPI_DOUBLE,targ[i],targ[i],comm,reqs++);
    start+=nsn;
  }
  for(reqs=info->reqs,i=np*2;i;--i) MPI_Wait(reqs++,&status);
  sh_ind = info->sh_ind;
  switch(op) {
    case OP_ADD: for(i=0;i<np;++i) { int c,j; for(c=nshared[i];c;--c) {
        double *uu=u+n*(*sh_ind++);
        for(j=n;j;--j) *uu++ += *buf++;
      } } break;
    case OP_MUL: for(i=0;i<np;++i) { int c,j; for(c=nshared[i];c;--c) {
        double *uu=u+n*(*sh_ind++);
        for(j=n;j;--j) *uu++ *= *buf++;
      } } break;
    case OP_MIN: for(i=0;i<np;++i) { int c,j; for(c=nshared[i];c;--c) {
        double *uu=u+n*(*sh_ind++);
        for(j=n;j;--j) { if(*buf < *uu) *uu = *buf; ++u, ++buf; }
      } } break;
    case OP_MAX: for(i=0;i<np;++i) { int c,j; for(c=nshared[i];c;--c) {
        double *uu=u+n*(*sh_ind++);
        for(j=n;j;--j) { if(*buf > *uu) *uu = *buf; ++u, ++buf; }
      } } break;
  }
}

static void nonlocal_transfer_many(
               double **u, int n, int op,
               const nonlocal_transfer_info *info, MPI_Comm comm)
{
  MPI_Status status;
  int np = info->np, i;
  MPI_Request *reqs = info->reqs;
  int *targ = info->target;
  int *nshared = info->nshared;
  int *sh_ind = info->sh_ind;
  int id;
  double *buf = info->buf, *start;
  MPI_Comm_rank(comm,&id);
  for(i=0;i<np;++i) {
    int c, j, ns = nshared[i];
    start = buf;
    for(j=0;j<n;++j) {double*uu=u[j]; for(c=0;c<ns;++c) *buf++=uu[sh_ind[c]];}
    sh_ind+=ns;
    MPI_Isend(start,n*ns,MPI_DOUBLE,targ[i],id,comm,reqs++);
  }
  start = buf;
  for(i=0;i<np;++i) {
    int nsn = n*nshared[i];
    MPI_Irecv(start,nsn,MPI_DOUBLE,targ[i],targ[i],comm,reqs++);
    start+=nsn;
  }
  for(reqs=info->reqs,i=np*2;i;--i) MPI_Wait(reqs++,&status);
  sh_ind = info->sh_ind;
  switch(op) {
    case OP_ADD: for(i=0;i<np;++i) { int c, j, ns = nshared[i];
      for(j=0;j<n;++j) { double *uu=u[j]; for(c=0;c<ns;++c)
        uu[sh_ind[c]] += *buf++;
      } sh_ind += ns; } break;
    case OP_MUL: for(i=0;i<np;++i) { int c, j, ns = nshared[i];
      for(j=0;j<n;++j) { double *uu=u[j]; for(c=0;c<ns;++c)
        uu[sh_ind[c]] *= *buf++;
      } sh_ind += ns; } break;
    case OP_MIN: for(i=0;i<np;++i) { int c, j, ns = nshared[i];
      for(j=0;j<n;++j) { double *uu=u[j]; for(c=0;c<ns;++c) {
        if(*buf < uu[sh_ind[c]]) uu[sh_ind[c]] = *buf;
        ++buf;
      } } sh_ind += ns; } break;
    case OP_MAX: for(i=0;i<np;++i) { int c, j, ns = nshared[i];
      for(j=0;j<n;++j) { double *uu=u[j]; for(c=0;c<ns;++c) {
        if(*buf > uu[sh_ind[c]]) uu[sh_ind[c]] = *buf;
        ++buf;
      } } sh_ind += ns; } break;
  }
}
#endif

/*--------------------------------------------------------------------------
   Combined Execution
  --------------------------------------------------------------------------*/

typedef struct {
  int *local_cm; /* local condense map */
#ifdef MPI
  nonlocal_transfer_info *nlinfo;
  MPI_Comm comm;
#endif
} transfer_info;

static void full_exec(double *u, int op, const transfer_info *info)
{
  local_condense(u,op,info->local_cm);
#ifdef MPI
  nonlocal_transfer(u,op,info->nlinfo,info->comm);
#endif
  local_uncondense(u,info->local_cm);
}

static void full_exec_vec(double *u, int n, int op, const transfer_info *info)
{
#ifdef MPI
  if(n>info->nlinfo->maxv)
    fail("gs2.c: initialized with max vec size = %d,"
         " but called with vec size = %d\n",info->nlinfo->maxv,n);
#endif
  local_condense_vec(u,n,op,info->local_cm);
#ifdef MPI
  nonlocal_transfer_vec(u,n,op,info->nlinfo,info->comm);
#endif
  local_uncondense_vec(u,n,info->local_cm);
}

static void full_exec_many(double **u, int n, int op, const transfer_info *info)
{
  int i;
#ifdef MPI
  if(n>info->nlinfo->maxv)
    fail("gs2.c: initialized with max vec size = %d,"
         " but called with vec size = %d\n",info->nlinfo->maxv,n);
#endif
  for(i=0;i<n;++i) local_condense(u[i],op,info->local_cm);
#ifdef MPI
  nonlocal_transfer_many(u,n,op,info->nlinfo,info->comm);
#endif
  for(i=0;i<n;++i) local_uncondense(u[i],info->local_cm);
}

/*--------------------------------------------------------------------------
   Crystal Router
  --------------------------------------------------------------------------*/

#ifdef MPI

typedef struct {
  int* data;
  int len;
  int max;
} buffer_type;

static void make_buffer(buffer_type *buffer, int max)
{
  buffer->data = (int*) smalloc(max*sizeof(int));
  buffer->max = max;
  buffer->len = 0;
}

/* split messages into those targeted for
       lo: dest id <  cutoff
       hi: dest id >= cutoff */
static void crystal_router_split(const buffer_type *src, int cutoff,
                                 buffer_type *lo, buffer_type *hi)
{
  int *buf_src = src->data;
  int buflen = src->len;

  lo->len = hi->len = 0;
  while(buflen>0) {
    int chlen = 3+buf_src[2];
    if(buf_src[0]<cutoff) {
      memcpy(lo->data+lo->len,buf_src,chlen*sizeof(int));
      lo->len+=chlen;
    } else {
      memcpy(hi->data+hi->len,buf_src,chlen*sizeof(int));
      hi->len+=chlen;
    }
    buf_src+=chlen, buflen-=chlen;
  }
}

/*
  exchange a sequence of (variable length integer sequence) messages
  with other processors, using the crystal routing algorithm
  
  the msgbuf argument points to a sequence of messages for both input and
  output, where each message has the following format:
  
     destination processor : 1 int
     source processor      : 1 int
     message length        : 1 int
     message               : <message length> ints
  
  arguments:
     buffer->data: pointer to message sequence; on input all messages should
                   have this processor as source; on output all messages will
                   have this processor as destination
     buffer->len:  total number of ints in message sequence on input and output
     buffer->max:  maximum buffer length (for output and intermediate stages)
*/
static void crystal_router(buffer_type *buffer, MPI_Comm comm)
{
  int base=0,num,id;
  buffer_type buf_lo, buf_hi;
  
  MPI_Comm_rank(comm,&id);
  MPI_Comm_size(comm,&num);
  if(num==1) return;

  make_buffer(&buf_lo,buffer->max);
  make_buffer(&buf_hi,buffer->max);

  while(num>1) {
    MPI_Request req1,req2,req3;
    MPI_Status status;
    int count, count2;
    
    /* we must communicate yet with procs [base,base+num) */
    int sz1 = num/2, cutoff=base+sz1;
    
    /* split messages into those targeted for
       lo: [base,  cutoff=base+sz1)
       hi: [cutoff,base+num) */
    crystal_router_split(buffer, cutoff, &buf_lo, &buf_hi);
    
    if(id<cutoff) {
      /* proc in lo group;
         with partner in hi group, trade our messages for hi group
                                   for their messages for lo group */
      int targ = id+sz1;  /* partner proc in hi group */
      
      /* if there is one more proc in the hi group than the lo group
         (i.e., num is odd), then the last proc of the lo group
         will receive messages from targ and targ+1 (the last proc of
         the hi group that would otherwise be partnerless) */
      int recvtwo = num&1 && id==base+sz1-1;
      
#ifdef DYNAMIC
      MPI_Isend(&buf_hi.len,1,MPI_INT,targ,  id,    comm,&req1);
      MPI_Irecv(&count,     1,MPI_INT,targ,  targ,  comm,&req2);
      if(recvtwo)
      MPI_Irecv(&count2,    1,MPI_INT,targ+1,targ+1,comm,&req3);
      MPI_Wait(&req1,&status); MPI_Wait(&req2,&status);
      if(recvtwo) { MPI_Wait(&req3,&status); count+=count2; }
      if(buf_lo.len+count > buf_lo.max) {
        count2 = buffer->max ? buffer->max : 1;
        while(count2 < buf_lo.len+count) count2*=2;
        free(buffer->data); make_buffer(buffer, count2);
        buf_lo.data = (int*)srealloc(buf_lo.data,count2*sizeof(int));
        buf_lo.max = count2;
        buf_hi.data = (int*)srealloc(buf_hi.data,count2*sizeof(int));
        buf_hi.max = count2;
      }
#endif

      MPI_Isend(buf_hi.data,buf_hi.len,MPI_INT,targ,id,comm,&req1);
      MPI_Irecv(buf_lo.data+buf_lo.len,buf_lo.max-buf_lo.len,MPI_INT,
                    targ,targ,comm,&req2);
      if(recvtwo) ++targ, MPI_Irecv(buffer->data,buffer->max,MPI_INT,
                                    targ,targ,comm,&req3);
      MPI_Wait(&req1,&status);
      MPI_Wait(&req2,&status);
      MPI_Get_count(&status,MPI_INT,&count);
      buf_lo.len+=count;
      if(recvtwo) {
        MPI_Wait(&req3,&status);
        MPI_Get_count(&status,MPI_INT,&count);
        if(buf_lo.len+count>buf_lo.max)
          failwith("gs2.c: insufficient buffer size in crystal_router");
        memcpy(buf_lo.data+buf_lo.len,buffer->data,count*sizeof(int));
        buf_lo.len+=count;
      }
      /* recurse */
      buffer->len=buf_lo.len;
      memcpy(buffer->data,buf_lo.data,buf_lo.len*sizeof(int));
      num=sz1;
    } else {
      /* proc in hi group;
         with partner in lo group, trade our messages for lo group
                                   for their messages for hi group */
      int targ = id-sz1;  /* partner proc in lo group */
       /* if targ is not in the lo group (because num was odd),
          then send our messages to the last proc of the lo group;
          additionally, we won't be receiving anything
       */
      int recving = 1;
      if(targ == cutoff) --targ, recving=0;

#ifdef DYNAMIC
      MPI_Isend(&buf_lo.len,1,MPI_INT,targ,  id,    comm,&req1);
      if(recving)
      MPI_Irecv(&count,     1,MPI_INT,targ,  targ,  comm,&req2);
      MPI_Wait(&req1,&status);
      if(recving) MPI_Wait(&req2,&status);
      if(recving && buf_hi.len+count > buf_hi.max) {
        count2 = buffer->max ? buffer->max : 1;
        while(count2 < buf_hi.len+count) count2*=2;
        free(buffer->data); make_buffer(buffer, count2);
        buf_lo.data = (int*)srealloc(buf_lo.data,count2*sizeof(int));
        buf_lo.max = count2;
        buf_hi.data = (int*)srealloc(buf_hi.data,count2*sizeof(int));
        buf_hi.max = count2;
      }
#endif

      MPI_Isend(buf_lo.data,buf_lo.len,MPI_INT,targ,id,comm,&req1);
      if(recving) MPI_Irecv(buf_hi.data+buf_hi.len,buf_hi.max-buf_hi.len,
                            MPI_INT,targ,targ,comm,&req2);
      MPI_Wait(&req1,&status);
      if(recving) {
        MPI_Wait(&req2,&status);
        MPI_Get_count(&status,MPI_INT,&count);
        buf_hi.len+=count;
      }
      /* recurse */
      buffer->len=buf_hi.len;
      memcpy(buffer->data,buf_hi.data,buf_hi.len*sizeof(int));
      base+=sz1, num-=sz1;
    }
  }
  free(buf_lo.data); free(buf_hi.data);
}
#endif

/*--------------------------------------------------------------------------
   Integer Map

   the type imap is a red-black tree taking an integer key to a linked list
   of integer values, which may be traversed in ascending order of keys
   
   [macro] void imap_new(imap *x, int size)
      allocates all memory needed for new map x,
      imap_insert(x, ...) may be called at most size times (no check)
   [macro] void imap_free(imap *x)
      deallocate memory associated with map x
   
   [macro] void imap_begin(imap_iterator p, imap *x)
      set p to the first node in the map
   [macro] void imap_next(imap_iterator p)
      advance p to the next node or to 0 (null) if p was the last node
      
   [macro] int imap_key(imap_iterator p)
      the key of node p
   [macro] void imap_ll_begin(imap_ll_iterator q, imap_iterator p)
      set q to the first linked list value node of map node p
   [macro] imap_ll_iterator imap_ll_next(imap_ll_iterator q)
      q->next, the next linked list value node (or 0)
   [macro] int imap_ll_val(imap_ll_iterator q)
      q->val, the value of node q

   [function] void imap_insert(imap *m, int key, int val)
      add (key, val) to the map
   [function] int imap_find_first(imap *m, int key)
      returns the most recent value associated with given key

  --------------------------------------------------------------------------*/

typedef struct imap_ll_node_ {
  struct imap_ll_node_ *next;
  int val;
} imap_ll_node;

typedef struct imap_rb_tree_node_ {
  struct imap_rb_tree_node_ *left, *right, *parent;
  int color;
  int key;
  imap_ll_node *vals;
} imap_rb_tree_node;

#define imap_color_black 0
#define imap_color_red 1

typedef struct {
  imap_rb_tree_node *root, *memt, *memtb;
  imap_ll_node *meml, *memlb;
} imap;

typedef imap_rb_tree_node *imap_iterator;
#define imap_begin(it,mp) { it=mp->root; if(it) while(it->left) it=it->left; }
#define imap_next(it) \
  { if(it->right) { it=it->right; while(it->left) it=it->left; } \
    else { imap_iterator oldimapiter=it; \
           while((it=it->parent)) { if(it->left==oldimapiter) break; \
                                  oldimapiter=it; } } }
#define imap_new(x,size) { x = (imap*)smalloc(sizeof(imap)); \
  x->root = 0, x->memt = x->memtb = (imap_rb_tree_node*) \
                        smalloc((size)*sizeof(imap_rb_tree_node)), \
               x->meml = x->memlb = (imap_ll_node*) \
                        smalloc((size)*sizeof(imap_ll_node)); }
#define imap_free(x) ( free(x->memlb), free(x->memtb), free(x) )
#define imap_key(it) ( (it)->key )

typedef imap_ll_node *imap_ll_iterator;
#define imap_ll_begin(it,mit) { it = mit->vals; }
#define imap_ll_next(it) { it=it->next; }
#define imap_ll_val(it) ( it->val )

/*
void imap_print_node(int lvl, imap_rb_tree_node* n)
{
  int i;
  imap_ll_iterator p;
  if(!n) return;
  for(i=0;i<lvl;++i) printf("  ");
  printf(n->color==imap_color_red?"r":"b");
  printf(" %d ->",n->key);
  imap_ll_begin(p,n);
  for(;p;p=p->next) printf(" %d", imap_ll_val(p));
  printf("\n");
  imap_print_node(lvl+1,n->left);
  imap_print_node(lvl+1,n->right);
}
*/

/* if the key is already in the map, val is added on to the
   linked list of values for the key */
static void imap_insert(imap *m, int key, int val)
{
  imap_rb_tree_node *p, *x;
  
  x = m->memt++;
  x->parent=x->left=x->right=0, x->color=imap_color_red, x->key=key;
  x->vals = m->meml++;
  x->vals->val = val, x->vals->next = 0;
  
  p = m->root;
  if(!p) { x->color=imap_color_black; m->root=x; return; }

  for(;;) {
    if(p->key==x->key) {
      x->vals->next=p->vals, p->vals=x->vals, --m->memt; return;
    } else if(x->key < p->key) {
      if(!p->left) { p->left=x, x->parent=p; break; }
      p=p->left;
    } else {
      if(!p->right) { p->right=x, x->parent=p; break; }
      p=p->right;
    }
  }

  for(;;) { /* check validity of node x */
    imap_rb_tree_node *p = x->parent, *g, *u, *gg;
    if(!p) { x->color = imap_color_black; return; } /* x is root */
    /* x is not root */
    if(p->color==imap_color_black) return;
    /* parent is red (and must itself have a parent) */
    g = p->parent; /* grand parent */
    u = p==g->left ? g->right : g->left; /* uncle */
    if(u && u->color==imap_color_red) {
      p->color = u->color = imap_color_black;
      g->color = imap_color_red;
      x = g;
      continue; /* recurse, check grandparent */
    }
    /* uncle is black (could be null) */
    if(x==p->right && p==g->left) {
      /* rotate left at p */
      p->right = x->left;
      if(p->right) p->right->parent = p;
      x->left = p; p->parent = x;
      g->left = x; x->parent = g;
      x = p;
      p = x->parent;
    } else if(x==p->left && p==g->right) {
      /* rotate right at p */
      p->left = x->right;
      if(p->left) p->left->parent = p;
      x->right = p; p->parent = x;
      g->right = x; x->parent = g;
      x = p;
      p = x->parent;
    }
    /* now either x==p->left && p==g->left or
                  x==p->right && p==g->right */
    p->color = imap_color_black;
    g->color = imap_color_red;
    gg = g->parent; /* great grand parent (may be null) */
    if(x==p->left) { /* && p==g->left */
      /* rotate right at g */
      g->left = p->right;
      if(g->left) g->left->parent = g;
      p->right = g; g->parent = p;
    } else { /* x==p->right && p==g->right */
      /* rotate left at g */
      g->right = p->left;
      if(g->right) g->right->parent = g;
      p->left = g; g->parent = p;
    }
    /* last step of the above rotations */
    if(gg) {
      if(gg->left==g) gg->left = p; else gg->right = p;
    }
    p->parent = gg;
    if(!gg) m->root=p;
    return;
  }
}

/* returns the first value associated with a given key
   (the first value in the linked list) */
static int imap_find_first(const imap *m, int key)
{
  imap_rb_tree_node *x;
  for(x=m->root;x;) {
    if(x->key==key) {
      if(!x->vals)
        failwith("gs2.c: logic error in discovery phase of gather/scatter");
      else return x->vals->val;
    } else
      x = key < x->key ? x->left : x->right;
  }
  failwith("gs2.c: logic error in discovery phase of gather/scatter");
  return 0;
}

/*--------------------------------------------------------------------------
   Integer Pair Set

   very similar to integer maps, above, but has no values, just keys,
   and the keys are integer pairs, with order by first then second integer
   
   [macro] void iset_new(iset *x, int size)
      allocates all memory needed for new set x,
      iset_insert(x, ...) may be called at most size times (no check)
   [macro] void imap_free(iset *x)
      deallocate memory associated with set x
   
   [macro] void iset_begin(iset_iterator p, iset *x)
      set p to the first node in the set
   [macro] void iset_next(iset_iterator p)
      advance p to the next node or to 0 (null) if p was the last node
      
   [macro] int iset_key1(iset_iterator p)
      the first integer of the key of node p
   [macro] int iset_key2(iset_iterator p)
      the second integer of the key of node p

   [function] void iset_insert(iset *x, int key1, int key2)
      add (key1, key2) to the set

  --------------------------------------------------------------------------*/

typedef struct iset_rb_tree_node_ {
  struct iset_rb_tree_node_ *left, *right, *parent;
  int color;
  int key1, key2;
} iset_rb_tree_node;

#define iset_color_black imap_color_black
#define iset_color_red imap_color_red

typedef struct iset_ {
  iset_rb_tree_node *root, *mem, *memb;
} iset;

typedef iset_rb_tree_node *iset_iterator;
#define iset_begin(it,set) imap_begin(it,set)
#define iset_next(it) \
  { if(it->right) { it=it->right; while(it->left) it=it->left; } \
    else { iset_iterator oldisetiter=it; \
           while((it=it->parent)) { if(it->left==oldisetiter) break; \
                                  oldisetiter=it; } } }
#define iset_new(x,size) { x = (iset*)smalloc(sizeof(iset)); \
  x->root = 0, x->mem = x->memb = (iset_rb_tree_node*) \
                        smalloc((size)*sizeof(iset_rb_tree_node)); }
#define iset_free(x) ( free(x->memb), free(x) )
#define iset_key1(it) ( (it)->key1 )
#define iset_key2(it) ( (it)->key2 )

/*
void iset_print_node(int lvl, iset_rb_tree_node* n)
{
  int i;
  if(!n) return;
  for(i=0;i<lvl;++i) printf("  ");
  printf(n->color==imap_color_red?"r":"b");
  printf(" (%d,%d)\n",n->key1,n->key2);
  iset_print_node(lvl+1,n->left);
  iset_print_node(lvl+1,n->right);
}
*/

static void iset_insert(iset *s, int key1, int key2)
{
  iset_rb_tree_node *p, *x;
  
  x = s->mem++;
  x->parent=x->left=x->right=0, x->color=iset_color_red;
  x->key1=key1, x->key2=key2;
  
  p = s->root;
  if(!p) { x->color=iset_color_black; s->root=x; return; }

  for(;;) {
    int cmp = x->key1 - p->key1;
    cmp = cmp ? cmp : x->key2 - p->key2;
    if(!cmp) {
      --s->mem; return;
    } else if(cmp < 0) {
      if(!p->left) { p->left=x, x->parent=p; break; }
      p=p->left;
    } else {
      if(!p->right) { p->right=x, x->parent=p; break; }
      p=p->right;
    }
  }

  for(;;) { /* check validity of node x */
    iset_rb_tree_node *p = x->parent, *g, *u, *gg;
    if(!p) { x->color = iset_color_black; return; } /* x is root */
    /* x is not root */
    if(p->color==iset_color_black) return;
    /* parent is red (and must itself have a parent) */
    g = p->parent; /* grand parent */
    u = p==g->left ? g->right : g->left; /* uncle */
    if(u && u->color==iset_color_red) {
      p->color = u->color = iset_color_black;
      g->color = iset_color_red;
      x = g;
      continue; /* recurse, check grandparent */
    }
    /* uncle is black (could be null) */
    if(x==p->right && p==g->left) {
      /* rotate left at p */
      p->right = x->left;
      if(p->right) p->right->parent = p;
      x->left = p; p->parent = x;
      g->left = x; x->parent = g;
      x = p;
      p = x->parent;
    } else if(x==p->left && p==g->right) {
      /* rotate right at p */
      p->left = x->right;
      if(p->left) p->left->parent = p;
      x->right = p; p->parent = x;
      g->right = x; x->parent = g;
      x = p;
      p = x->parent;
    }
    /* now either x==p->left && p==g->left or
                  x==p->right && p==g->right */
    p->color = iset_color_black;
    g->color = iset_color_red;
    gg = g->parent; /* great grand parent (may be null) */
    if(x==p->left) { /* && p==g->left */
      /* rotate right at g */
      g->left = p->right;
      if(g->left) g->left->parent = g;
      p->right = g; g->parent = p;
    } else { /* x==p->right && p==g->right */
      /* rotate left at g */
      g->right = p->left;
      if(g->right) g->right->parent = g;
      p->left = g; g->parent = p;
    }
    /* last step of the above rotations */
    if(gg) {
      if(gg->left==g) gg->left = p; else gg->right = p;
    }
    p->parent = gg;
    if(!gg) s->root=p;
    return;
  }
}

/*--------------------------------------------------------------------------
   Computation of Transfer Information
  --------------------------------------------------------------------------*/

/* computes (and allocates) the local condense map given a map
   from global to local indices */
static int *mk_local_cm(const imap *m)
{
  int count = 1, *cm, *q;
  imap_iterator p;
  imap_begin(p,m);
  while(p) {
    imap_ll_iterator r;
    imap_ll_begin(r,p);
    if(r->next) {
      ++count;
      while(r) ++count, r=r->next;
    }
    imap_next(p);
  }
  q = cm = (int*)smalloc(count*sizeof(int));
  imap_begin(p,m);
  while(p) {
    imap_ll_iterator r;
    imap_ll_begin(r,p);
    if(r->next) {
      while(r) *q++ = r->val, r=r->next;
      *q++ = -1;
    }
    imap_next(p);
  }
  *q = -1;
  return cm;
}

#ifdef MPI   

/* given map from global to local indices (on this proc),
   receives a list of all points such that
     (global index) mod (number of processors) = (this proc id)
   employing the crystal router */
static void bin_from_glm(const imap *glm, buffer_type *buf, MPI_Comm comm)
{
  int buflen = 0;
  int *binsize, *binoffset, offset=0, *data;
  int id, np, i;
  
  imap_iterator p;
  
  MPI_Comm_rank(comm,&id);
  MPI_Comm_size(comm,&np);
  
  binsize = (int*)scalloc(2*np,sizeof(int));
  binoffset = binsize+np;
  
  imap_begin(p,glm);
  while(p) {
    ++binsize[imap_key(p)%np];
    imap_next(p);
  }
  
  for(i=0;i<np;++i) {
    if(!binsize[i]) continue;
    buflen+=3+binsize[i];
    binoffset[i]=buflen;
  }
  if(buflen>buf->max) {
    free(buf->data);
    make_buffer(buf,buflen);
  }
  buf->len=buflen, data=buf->data;
  for(i=0;i<np;++i) {
    if(!binsize[i]) continue;
    data[offset+0] = i; /* message destination */
    data[offset+1] = id; /* message source */
    data[offset+2] = binsize[i]; /* message length */
    offset+=3+binsize[i];
  }
  
  imap_begin(p,glm);
  while(p) {
    int gli = imap_key(p);
    data[--binoffset[gli%np]] = gli;
    imap_next(p);
  }
  
  free(binsize);
  crystal_router(buf,comm);
}

/* sorts the data from previous function by making a map of
   (global index) to (processor list) */
static imap *gpm_from_bin(const buffer_type *buffer)
{
  int len = buffer->len, *buf = buffer->data;
  imap *gpm;
  imap_new(gpm,2*len);
  while(len>0) {
    int src=buf[1], n = buf[2];
    buf+=3,len-=3+n;
    for(;n;--n) imap_insert(gpm, *buf++, src);
  }
  return gpm;
}

/* rearrange data (sorting) from previous function into a map of
     (processor 1, processor 2) to (global index list)
   where procs 1 & 2 share the list of global indices
   implemented by associating (proc1,proc2) with proc1*num_prcs+proc2 */
static imap *ppgm_from_gpm(const imap *gpm, int np)
{
  imap *ppgm;
  imap_iterator p;
  int count = 0;
  imap_begin(p,gpm);
  while(p) {
    int c = 0;
    imap_ll_iterator q;
    imap_ll_begin(q,p);
    for(;q;q=q->next) count+=c, c+=2;
    imap_next(p);
  }
  
  imap_new(ppgm,count);
  imap_begin(p,gpm);
  while(p) {
    int gli = imap_key(p);
    imap_ll_iterator q, r;
    imap_ll_begin(q,p);
    for(;q;q=q->next) {
      int p1 = imap_ll_val(q);
      for(r=q->next;r;r=r->next) {
        int p2 = imap_ll_val(r);
        imap_insert(ppgm,p1*np+p2,gli);
        imap_insert(ppgm,p2*np+p1,gli);
      }
    }
    imap_next(p);
  }
  
  return ppgm;
}

/* pack the data in the map from the above function
   into a sequence of messages for other procs
   exchange messages with the crystal router */
static void slbuf_from_ppgm(const imap *ppgm, buffer_type *buffer,
                            MPI_Comm comm)
{
  int buflen = 0;
  int *binsize, **bin, offset=0, *buf;
  int id, np, i;
  
  imap_iterator p;
  
  MPI_Comm_rank(comm,&id);
  MPI_Comm_size(comm,&np);
  
  binsize = (int*)scalloc(np,sizeof(int));
  bin = (int**)scalloc(np,sizeof(int*));
  
  imap_begin(p,ppgm);
  while(p) {
    int p1 = imap_key(p)/np;
    imap_ll_iterator q;
    binsize[p1]+=2;
    imap_ll_begin(q,p);
    for(;q;q=q->next) ++binsize[p1];
    imap_next(p);
  }
  
  for(i=0;i<np;++i) {
    if(!binsize[i]) continue;
    buflen+=3+binsize[i];
  }
  if(buflen>buffer->max) {
    free(buffer->data);
    make_buffer(buffer,buflen);
  }
  buf = buffer->data, buffer->len = buflen;
  for(i=0;i<np;++i) {
    if(!binsize[i]) continue;
    buf[offset+0] = i; /* message destination */
    buf[offset+1] = id; /* message source */
    buf[offset+2] = binsize[i]; /* message length */
    bin[i] = buffer->data+offset+3;
    offset+=3+binsize[i];
    binsize[i]=0;
  }
  
  imap_begin(p,ppgm);
  while(p) {
    int pp = imap_key(p), p1 = pp/np, p2 = pp%np;
    int *q = bin[p1], count=0;
    imap_ll_iterator r;
    *q++ = p2;
    imap_ll_begin(r,p);
    for(;r;r=r->next) ++count;
    *q++ = count;
    imap_ll_begin(r,p);
    for(;r;r=r->next) *q++ = imap_ll_val(r);
    bin[p1]=q;
    imap_next(p);
  }
  
  free(binsize); free(bin);
  crystal_router(buffer,comm);
}

/* unpack the messages from other procs into a set of
   (processor, global index) pairs, each indicating
   that we share a global index with some other proc */
static iset *pgs_from_slbuf(const buffer_type *buffer)
{
  int *buf = buffer->data, buflen = buffer->len;
  iset *s;
  iset_new(s,buflen*2);
  while(buflen>0) {
    int len = buf[2];
    buflen-=len+3, buf+=3;
    while(len>0) {
      int p = *buf++;
      int n = *buf++;
      len-=n+2;
      for(;n;--n) iset_insert(s,p,*buf++);
    }
  }
  return s;
}

/* use the set from the above function to build the structure
   used for the execution phase */
static nonlocal_transfer_info *nlinfo_from_pgs(const iset *pgs,
                                               const imap *glm, int maxv)
{
  nonlocal_transfer_info *info =
    (nonlocal_transfer_info*) smalloc(sizeof(nonlocal_transfer_info));
  int count=0, oldp=-1,i,*sh_ind,np=0;
  iset_iterator q;
  
  iset_begin(q,pgs);
  while(q) {
    int p = iset_key1(q);
    if(p!=oldp) ++np, oldp=p;
    ++count;
    iset_next(q);
  }
  info->np = np;
  info->target = (int*)smalloc((2*np+count)*sizeof(int));
  info->nshared = info->target + np;
  sh_ind = info->sh_ind = info->nshared + np;
  info->reqs = (MPI_Request*)smalloc(2*np*sizeof(MPI_Request));
  info->buflen = 2*maxv*count; /* count for send, count for recv */
  info->buf = (double*)smalloc(info->buflen*sizeof(double));
  info->maxv = maxv;

  oldp=-1,i=-1;
  iset_begin(q,pgs);
  while(q) {
    int p = iset_key1(q), gli = iset_key2(q);
    if(p!=oldp) {
      oldp=p;
      if(i>=0) info->nshared[i]=count;
      info->target[++i]=p;
      count=0;
    }
    *sh_ind++ = imap_find_first(glm,gli);
    ++count;
    iset_next(q);
  }
  if(i>=0) info->nshared[i]=count;
  
  return info;
}

static void nlinfo_free(nonlocal_transfer_info *info)
{
  free(info->buf);
  free(info->reqs);
  free(info->target);
  free(info);
}

static nonlocal_transfer_info *nlinfo_from_glm(const imap *glm, int maxbuf,
                                               int maxv, MPI_Comm comm)
{
  buffer_type buffer;
  int id,np;
  imap *gpm, *ppgm;
  iset *pgs;
  nonlocal_transfer_info *info;
  
  make_buffer(&buffer,maxbuf);
  MPI_Comm_size(comm,&np);
  MPI_Comm_rank(comm,&id);
  /* glm, given, is a map from
       global vertex index -> local vertex indices */
  bin_from_glm(glm,&buffer,comm);
  gpm = gpm_from_bin(&buffer);
  /* gpm is a map from
       global vertex index -> processors w/ this vertex
     where (global vertex index) mod np = id */
  ppgm = ppgm_from_gpm(gpm,np);
  imap_free(gpm);
  /* ppgm is a map from
       proc1 * np + proc2  -> vertices shared by proc1 and proc2
     the map is symmetric: entries for (proc1,proc2) and (proc2,proc1) are equal
     vertices listed are all congruent to id mod np */
  slbuf_from_ppgm(ppgm,&buffer,comm);
  imap_free(ppgm);
  pgs = pgs_from_slbuf(&buffer);
  /* pgs is a set of
       (processor number, global index)
     listing which global vertices we share with which processors */
  info = nlinfo_from_pgs(pgs,glm,maxv);
  iset_free(pgs);
  free(buffer.data);
  return info;
}

#endif

static void tinfo_from_va(const int *v, int vn, int maxv, transfer_info *info)
{
  imap *glm;
  int count=0,n,num;
  const int *p;
  for(p=v,n=vn;n;--n) if(*p++) ++count; /* count non-zero entries */
  
  imap_new(glm,count);
  for(n=0;n<vn;++n) if(v[n]) imap_insert(glm,v[n],n);
  
  info->local_cm = mk_local_cm(glm);
  /* notice that INITIAL_BUFFER_SIZE is chosen (somewhat arbitrarily) as
     the maximum buffer size below */
#ifdef MPI
  MPI_Comm_size(info->comm,&num);
  info->nlinfo = nlinfo_from_glm(glm,INITIAL_BUFFER_SIZE,maxv,info->comm);
#endif
  imap_free(glm);
}

/*--------------------------------------------------------------------------
   FORTRAN Interface
  --------------------------------------------------------------------------*/

static transfer_info *cpgs_info = 0;
static int cpgs_max = 0;
static int cpgs_n = 0;

#ifndef MPI
typedef void MPI_Comm;
#endif
void FORTRAN_SETUP_NAME(const int *v, const int *vn, const int *maxv,
                        const MPI_Comm *comm, const int *np, int *handle)
{
  int mv = *maxv <= 1 ? 1 : *maxv;
#ifdef MPI
  int mpi_np;
  MPI_Comm_size(*comm,&mpi_np);
  if(mpi_np!=*np)
    fail("cpgs_setup: passed np of %d, but MPI_Comm_size gives np=%d\n",
         *np, mpi_np);
#else
  if(*np!=1)
    fail("cpgs_setup: not compiled with -DMPI; can't handle np=%d\n", *np);
#endif
  if(cpgs_max==cpgs_n) {
    cpgs_max = cpgs_max?cpgs_max+(cpgs_max+1)/2:4;
    cpgs_info = (transfer_info*)
      srealloc(cpgs_info, cpgs_max*sizeof(transfer_info));
  }
#ifdef MPI
  MPI_Comm_dup(*comm,&cpgs_info[cpgs_n].comm);
#endif
  tinfo_from_va(v,*vn,mv,cpgs_info+cpgs_n);
  *handle = cpgs_n++;
}

void FORTRAN_OP_NAME(const int *handle, double *u, const int *op)
{
  if(*op<1 || *op>4) failwith("invalid operation to cgps_op");
  if(*handle<0 || *handle>=cpgs_n || !cpgs_info[*handle].local_cm)
    failwith("invalid handle to cgps_op");
  full_exec(u,*op,cpgs_info + *handle);
}

void FORTRAN_OP_VEC_NAME(const int *handle, double *u, const int *n,
                         const int *op)
{
  if(*op<1 || *op>4) failwith("invalid operation to cgps_op");
  if(*handle<0 || *handle>=cpgs_n || !cpgs_info[*handle].local_cm)
    failwith("invalid handle to cgps_op_vec");
  full_exec_vec(u,*n,*op,cpgs_info + *handle);
}

void FORTRAN_OP_MANY_NAME(const int *handle, double *u, double *v, double *w,
                          const int *n, const int *op)
{
  double *uu[3];
  uu[0]=u,uu[1]=v,uu[2]=w;
  if(*op<1 || *op>4) failwith("invalid operation to cgps_op");
  if(*handle<0 || *handle>=cpgs_n || !cpgs_info[*handle].local_cm)
    failwith("invalid handle to cgps_op_many");
  full_exec_many(uu,*n,*op,cpgs_info + *handle);
}

/* note that we let some bytes linger in the cpgs_info array ...
   can't rearrange the array, because the handles are just indices
   and would become invalid
   using a handle after it has been freed will result in undefined
   behavior
*/
void FORTRAN_FREE_NAME(int *handle)
{
  if(*handle<0 || *handle>=cpgs_n || !cpgs_info[*handle].local_cm)
    failwith("invalid handle to cgps_op");
  free(cpgs_info[*handle].local_cm);
#ifdef MPI
  nlinfo_free(cpgs_info[*handle].nlinfo);
  MPI_Comm_free(&cpgs_info[*handle].comm);
#endif
  cpgs_info[*handle].local_cm = 0;
}
