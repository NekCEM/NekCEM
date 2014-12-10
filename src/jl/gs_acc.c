#ifdef _OPENACC

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "c99.h"
#include "name.h"
#include "fail.h"
#include "types.h"
#define gs_op gs_op_t
#include "gs_defs.h"
#include "gs_local.h"
#include "comm.h"
#include "mem.h"
#include "sort.h"
#include "crystal.h"
#include "sarray_sort.h"
#include "sarray_transfer.h"
#define GS_ACC_C
#include "gs_acc.h"

GS_DEFINE_DOM_SIZES()

typedef enum { mode_plain, mode_vec, mode_many, mode_dry_run } gs_mode;

static buffer static_buffer = null_buffer;

struct pw_comm_data {
  uint n;      /* number of messages */
  uint *p;     /* message source/dest proc */
  uint *size;  /* size of message */
  uint total;  /* sum of message sizes */
};

struct pw_data {
  struct pw_comm_data comm[2];
  const uint *map[2];
  comm_req *req;
  uint buffer_size;
};

typedef void exec_fun(void *data, gs_mode mode, unsigned vn, gs_dom dom, gs_op op,
		      unsigned transpose, const void *execdata, const struct comm *comm, char *buf);

typedef void fin_fun(void *data);

struct gs_remote {
  uint buffer_size, mem_size;
  void *data;
  exec_fun *exec;
  fin_fun *fin;
};

struct gs_data {
  struct comm comm;
  const uint *map_local[2]; /* 0=unflagged, 1=all */
  const uint *flagged_primaries;
  struct gs_remote r;
  uint handle_size;
};

static char *pw_exec_recvs(char *buf, const unsigned unit_size, const struct comm *comm,
                           const struct pw_comm_data *c, comm_req *req)
{
  const uint *p, *pe, *size=c->size;
  for(p=c->p,pe=p+c->n;p!=pe;++p) {
    size_t len = *(size++)*unit_size;
    comm_irecv(req++,comm,buf,len,*p,*p);
    buf += len;
  }
  return buf;
}

static char *pw_exec_sends(char *buf, const unsigned unit_size, const struct comm *comm,
                           const struct pw_comm_data *c, comm_req *req)
{
  const uint *p, *pe, *size=c->size;
  for(p=c->p,pe=p+c->n;p!=pe;++p) {
    size_t len = *(size++)*unit_size;
    comm_isend(req++,comm,buf,len,*p,comm->id);
    buf += len;
  }
  return buf;
}

//
// The above is duplicated from gs.c
// The below is our "new" code
//

#include <openacc.h>

#define USE_GPU_DIRECT 1

static char *pw_exec_recvs_acc(char *buf, const unsigned unit_size, const struct comm *comm,
			       const struct pw_comm_data *c, comm_req *req, uint *nr)
{
  const uint *p, *pe, *size=c->size;
  uint  l=0;

  // Find the size used for buf so we can make it present()
  for(p=c->p,pe=p+c->n;p!=pe;++p) {
    l += *(size++)*unit_size;
  }
  size=c->size;

  // Now that we have the size, send with GPU-Direct
#pragma acc data present(buf[0:l])
  {
    l = 0;
    for(p=c->p,pe=p+c->n;p!=pe;++p) {
      size_t len = *(size++)*unit_size;
#pragma acc host_data use_device(buf)
      {
	if(len) {
	  comm_irecv(req++,comm,&(buf[l]),len,*p,*p);
	  (*nr)++;
	}
      }
      l += len;
    }
  }
  return buf;
}

static char *pw_exec_sends_acc(char *buf, const unsigned unit_size, const struct comm *comm,
			       const struct pw_comm_data *c, comm_req *req, uint *nr)
{
  const uint *p, *pe, *size=c->size;
  uint l=0;

  // Find the size used for buf so we can make it present()
  for(p=c->p,pe=p+c->n;p!=pe;++p) {
    l += *(size++)*unit_size;
  }
  size=c->size;

  // Now that we have the size, send with GPU-Direct
#pragma acc data present(buf[0:l])
  {
    l = 0;
    for(p=c->p,pe=p+c->n;p!=pe;++p) {
      size_t len = *(size++)*unit_size;
#pragma acc host_data use_device(buf)
      {
	if(len) {
	  comm_isend(req++,comm,&(buf[l]),len,*p,comm->id);
	  (*nr)++;
	}
      }
      l += len;
    }
  }
  return buf;
}

static int map_size(int *map, int *t)
{
  int i,ct=0;

  *t = 0;

  // No map
  if(!map) {
    return 0;
  }
  
  // "Empty" map (contains only a single -1 terminator)
  if(map[0] == -1) {
    return 1;
  }

  // "Regular" map (contains two -1's as termination)
  for(i=ct=0;ct<2;i++){
    if(map[i]==-1){
      ct++;
      (*t)++;
    } else {
      ct=0;
    }
  }
  (*t)--;
  printf("");
  return i;
}

void fgs_fields_acc(const sint *handle, double *u, const sint *stride, const sint *n,
		    const sint *dom, const sint *op, const sint *transpose,
		    struct gs_data **fgs_info)
{
  struct pw_data *pwd;
  struct comm    *comm;
  buffer         *buf;
  const unsigned recv = 0^*transpose, send = 1^*transpose;
  uint    i,j,k;
  uint    bs,bl,uds,dstride,dtrans,vn,nr;
  uint    m_size,fp_m_size,snd_m_size,rcv_m_size,t_m_size;
  int     m_nt,fp_m_nt,snd_m_nt,rcv_m_nt,t_m_nt;
  int    *map,*t_map,*fp_map,*snd_map,*rcv_map;
  int    *mapf,*t_mapf,*fp_mapf,*snd_mapf,*rcv_mapf;
  double  *dbufp,*sbuf,*rbuf;
  double  t;

  char   hname[1024];
  static int calls=0;

  // Flatten...
  dstride = *stride;
  dtrans  = *transpose;
  vn      = *n;

  // Create temp buffer for gather/scatter and send/recv
  buf = &static_buffer;
  bs = (*n)*sizeof(double)*(fgs_info[*handle]->r.buffer_size);
  bl = (bs/sizeof(double))/2;
  //buffer_reserve(buf,bs);
  //dbufp   = (double*)buf->ptr;
  //sbuf    = dbufp;
  //rbuf    = &(dbufp[bl]);
  sbuf    = malloc(bl*sizeof(double));
  rbuf    = malloc(bl*sizeof(double));
  uds     = (*stride) * (*n); // Get size of u in number of doubles
  pwd     = fgs_info[*handle]->r.data;
  comm    = &fgs_info[*handle]->comm;

  // Flatten...
  map        = (int*)(fgs_info[*handle]->map_local[0^*transpose]);
  t_map      = (int*)(fgs_info[*handle]->map_local[1^*transpose]);
  fp_map     = (int*)(fgs_info[*handle]->flagged_primaries);
  snd_map    = (int*)(pwd->map[send]);
  rcv_map    = (int*)(pwd->map[recv]);
  fp_m_size  = map_size(fp_map,&fp_m_nt);
  m_size     = map_size(map,&m_nt);  
  snd_m_size = map_size(snd_map,&snd_m_nt);
  rcv_m_size = map_size(rcv_map,&rcv_m_nt);
  t_m_size   = map_size(t_map,&t_m_nt);

  //If *stride is 0, we are getting called from gs_op
  // As such, we need to give a uds a value of 1 - Matt Otten 10/31
  if((*stride)==0) uds        = 1;

  mapf = (int*)malloc(m_nt*2*sizeof(int));
  for(i=0,k=0;map[i]!=-1;i=j+1,k++){
      // Recortd i
    mapf[k*2] = i;
      for(j=i+1;map[j]!=-1;j++);
      // Record j-i
      mapf[k*2+1] = j-i-1;
  }

  t_mapf = (int*)malloc(t_m_nt*2*sizeof(int));
  for(i=0,k=0;t_map[i]!=-1;i=j+1,k++){
      // Recortd i
      t_mapf[k*2] = i;
      for(j=i+1;t_map[j]!=-1;j++);
      // Record j-i
      t_mapf[k*2+1] = j-i-1;
  }

  fp_mapf = (int*)malloc(fp_m_nt*2*sizeof(int));
  for(k=0;k<fp_m_nt;k++){
    // Recortd i
    fp_mapf[k*2] = 0;
    for(i=0;fp_map[i]!=-1;i++);
    // Record j-i
    fp_mapf[k*2+1] = i;
  }

  snd_mapf = (int*)malloc(snd_m_nt*2*sizeof(int));
  for(i=0,k=0;snd_map[i]!=-1;i=j+1,k++){
    // Recortd i
    snd_mapf[k*2] = i;
    for(j=i+1;snd_map[j]!=-1;j++);
    // Record j-i
    snd_mapf[k*2+1] = j-i-1;
  }

  rcv_mapf = (int*)malloc(rcv_m_nt*2*sizeof(int));
  for(i=0,k=0;rcv_map[i]!=-1;i=j+1,k++){
    // Recortd i
    rcv_mapf[k*2] = i;
    for(j=i+1;rcv_map[j]!=-1;j++);
    // Record j-i
    rcv_mapf[k*2+1] = j-i-1;
  }

#if 1
  calls++;
  gethostname(hname, sizeof(hname));
  //fprintf(stderr,"%s: enter %d\n",hname,calls);
#endif
#if 0
#endif
#if 0
  fprintf(stderr,"%s: map[0:%d]     -> %lX : %lX\n",hname,m_size,map,map+m_size);
  fprintf(stderr,"%s: t_map[0:%d]   -> %lX : %lX\n",hname,t_m_size,t_map,t_map+t_m_size);
  fprintf(stderr,"%s: fp_map[0:%d]  -> %lX : %lX\n",hname,fp_m_size,fp_map,fp_map+fp_m_size);
  fprintf(stderr,"%s: snd_map[0:%d] -> %lX : %lX\n",hname,snd_m_size,snd_map,snd_map+snd_m_size);
  fprintf(stderr,"%s: rcv_map[0:%d] -> %lX : %lX\n",hname,rcv_m_size,rcv_map,rcv_map+rcv_m_size);
#endif

#pragma acc enter data pcopyin(t_map[0:t_m_size],map[0:m_size],fp_map[0:fp_m_size],snd_map[0:snd_m_size],rcv_map[0:rcv_m_size])
#pragma acc data present(u[0:uds]) copyin(t_mapf[0:t_m_nt*2],mapf[0:m_nt*2],snd_mapf[0:snd_m_nt*2],rcv_mapf[0:rcv_m_nt*2])
  {
#pragma acc data create(sbuf[0:bl],rbuf[0:bl]) if(bl!=0)
    {
      // The below implementing cgs_many()/gs_aux():
      //
      // gs_aux_acc(u,mode_many,dn,fgs_dom[*dom],(gs_op_t)(*op-1),*transpose!=0,fgs_info[*handle],NULL,us);
      //
      
      // gs_gather_many_acc(u,u,vn,gsh->map_local[0^transpose],dom,op); 
      for(k=0;k<vn;++k) {
#pragma acc parallel loop gang vector present(u[0:uds],map[0:m_size],mapf[0:m_nt*2]) private(i,j,t) async(k)
	for(i=0;i<m_nt;i++){
	  t = u[map[mapf[i*2]]+k*dstride];
#pragma acc loop seq
	  for(j=0;j<mapf[i*2+1];j++) {
	    t += u[map[mapf[i*2]+j+1]+k*dstride];
	  }
	  u[map[mapf[i*2]]+k*dstride] = t;
	}
      }
#pragma acc wait      
      /*
      for(k=0;k<vn;++k) {
	for(i=0;map[i]!=-1;i=j+1){
	  t = u[map[i]+k*dstride];
	  for(j=i+1;map[j]!=-1;j++){
	    t += u[map[j]+k*dstride];
	  }
	  u[map[i]+k*dstride] = t;
	}
      }
      */
      
      if(dtrans==0) {
        // gs_init_many_acc(u,vn,gsh->flagged_primaries,dom,op);
        for(k=0;k<vn;++k) {
#pragma acc parallel loop gang vector present(u[0:uds],fp_map[0:fp_m_size],fp_mapf[0:fp_m_nt*2]) private(i,j) async(k)
          for(i=0;i<fp_m_nt;i++){
#pragma acc loop seq
            for(j=0;j<fp_mapf[i*2+1];j++){
              u[fp_map[fp_mapf[i*2]+j]+k*dstride]=0.0;
            }
          }
        }
#pragma acc wait
      }

      // --
      /*      if(dtrans==0) {
	// gs_init_many_acc(u,vn,gsh->flagged_primaries,dom,op);
#pragma acc parallel loop gang vector present(u[0:uds],fp_map[0:fp_m_size]) private(i,k)
	for(k=0;k<vn;++k) {
	  for(i=0;fp_map[i]!=-1;i++){
	    u[fp_map[i]+k*dstride]=0.0;
	  }
	}
      }
      */
      /* post receives */
      // If the send/recv buffer is of length 0, we skip all the MPI and asssociated pack/unpack.
      if( bl ) { 
#if USE_GPU_DIRECT
	nr = 0;
	pw_exec_recvs_acc((char*)rbuf,vn*sizeof(double),comm,&pwd->comm[recv],pwd->req,&nr);
#else
	pw_exec_recvs((char*)rbuf,vn*sizeof(double),comm,&pwd->comm[recv],pwd->req);
#endif
	/* fill send buffer */
	// gs_scatter_many_to_vec_acc(sendbuf,data,vn,pwd->map[send],dom);
#pragma acc parallel loop gang vector present(u[0:uds],snd_map[0:snd_m_size],snd_mapf[0:snd_m_nt*2],sbuf[0:bl]) private(i,j,k,t) collapse(2)
	for(k=0;k<vn;++k) {
	  for(i=0;i<snd_m_nt;i++){
#pragma acc loop seq
	    for(j=0;j<snd_mapf[i*2+1];j++) {
	      sbuf[k+snd_map[snd_mapf[i*2]+j+1]*vn] = u[snd_map[snd_mapf[i*2]]+k*dstride];
	    }
	  }
	}
	//#pragma acc wait      
	/*
#pragma acc parallel loop gang vector present(u[0:uds],snd_map[0:snd_m_size],sbuf[0:bl]) private(i,j,k)
	for(k=0;k<vn;k++) {
	  for(i=0;snd_map[i]!=-1;i=j+1){
	    for(j=i+1;snd_map[j]!=-1;j++){
	      //	      sbuf[k*dstride+snd_map[j]] = u[snd_map[i]+k*dstride];
	      sbuf[k+snd_map[j]*vn] = u[snd_map[i]+k*dstride];
	    }
	  }
	}
	*/
	/* post sends */
#if USE_GPU_DIRECT
	pw_exec_sends_acc((char*)sbuf,vn*sizeof(double),comm,&pwd->comm[send],&pwd->req[nr],&nr);
	comm_wait(pwd->req,nr);
#else
#pragma acc update host(sbuf[0:bl]) async(1)
#pragma acc wait	
	pw_exec_sends((char*)sbuf,vn*sizeof(double),comm,&pwd->comm[send],&pwd->req[pwd->comm[recv].n]);
	comm_wait(pwd->req,pwd->comm[0].n+pwd->comm[1].n);
#pragma acc update device(rbuf[0:bl]) async(1)
#pragma acc wait	
#endif
	/* gather using recv buffer */
	// gs_gather_vec_to_many_acc(data,buf,vn,pwd->map[recv],dom,op);
#pragma acc parallel loop gang vector present(u[0:uds],rcv_map[0:rcv_m_size],rcv_mapf[0:rcv_m_nt*2],rbuf[0:bl]) private(i,j,k,t) collapse(2)
	for(k=0;k<vn;++k) {
	  for(i=0;i<rcv_m_nt;i++){
#pragma acc loop seq
	    for(j=0;j<rcv_mapf[i*2+1];j++) {
	      u[rcv_map[rcv_mapf[i*2]]+k*dstride] += rbuf[k+rcv_map[rcv_mapf[i*2]+j+1]*vn];
	    }
	  }
	}
	//#pragma acc wait      
	/*
	#pragma acc parallel loop gang vector present(u[0:uds],rcv_map[0:rcv_m_size],rbuf[0:bl]) private(i,j,k)
	for(k=0;k<vn;k++){
	  for(i=0;rcv_map[i]!=-1;i=j+1){
	    for(j=i+1;rcv_map[j]!=-1;j++){
	      //	      u[rcv_map[i]+k*dstride] += rbuf[k*dstride+rcv_map[j]];
	      u[rcv_map[i]+k*dstride] += rbuf[k+rcv_map[j]*vn];
	    }
	  }
	}
	*/
      }
      // --
      // gs_scatter_many_acc(u,u,vn,gsh->map_local[1^transpose],dom); 
      for(k=0;k<vn;++k) {
#pragma acc parallel loop gang vector present(u[0:uds],t_map[0:t_m_size],t_mapf[0:t_m_nt*2]) private(i,j,t) async(k)
	for(i=0;i<t_m_nt;i++){
	  t = u[t_map[t_mapf[i*2]]+k*dstride];
#pragma acc loop seq
	  for(j=0;j<t_mapf[i*2+1];j++) {
	    u[t_map[t_mapf[i*2]+j+1]+k*dstride] = t;
	  }
	}
      }
#pragma acc wait      
      /*
      for(k=0;k<vn;++k) {
	for(i=0;t_map[i]!=-1;i=j+1) {
	  t=u[t_map[i]+k*dstride];
	  for(j=i+1;t_map[j]!=-1;j++) {
	    u[t_map[j]+k*dstride] = t;
	  }
	}
      }
      */
      // --
    }
  }
  // --
  if( bl ) {
    free(sbuf);
    free(rbuf);
  }
  free(mapf);
  free(t_mapf);
  free(fp_mapf);
  free(snd_mapf);
  free(rcv_mapf);
#if 0
  fprintf(stderr,"%s: exit %d\n",hname,calls);
#endif
}
#endif
