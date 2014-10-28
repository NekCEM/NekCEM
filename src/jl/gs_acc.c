#include <openacc.h>

void fgs_fields_acc(const sint *handle,
		    void *u, const sint *stride, const sint *n,
		    const sint *dom, const sint *op, const sint *transpose)
{
  struct pw_data *pwd;
  struct comm    *comm;
  buffer         *buf;
  const unsigned recv = 0^*transpose, send = 1^*transpose;
  uint    i,j,k,bs,uds,dstride,dtrans,vn;
  int    *l_map,*map;
  char   *bufp,*sbufp,*q;
  double  t,*ud;
  gs_dom  l_dom;

  fgs_check_parms(*handle,*dom,*op,"gs_op_fields",__LINE__);
  if(*n<0) return;

  // Setup a "double" version of u
  ud = u;                 // Get pointer from u
  uds = (*stride) * (*n); // Get size of u in number of doubles

  // Get the map and other pointers we need.
  map     = (int*)(fgs_info[*handle]->map_local[1^*transpose]);
  dstride = *stride;
  l_dom   = fgs_dom[*dom];
  pwd     = fgs_info[*handle]->r.data;
  comm    = &(fgs_info[*handle]->comm);
  dtrans  = *transpose;
  vn      = *n;

  // Create temp buffer for gather/scatter and send/recv
  buf = &static_buffer;
  bs = vn*gs_dom_size[l_dom]*(fgs_info[*handle]->r.buffer_size);
  buffer_reserve(buf,bs);
  bufp = buf->ptr;

  //#pragma acc data create(bufp[0:bs]) present(ud[0:uds])
  {
    // The below implementing cgs_many()/gs_aux():
    //
    // gs_aux_acc(u,mode_many,dn,fgs_dom[*dom],(gs_op_t)(*op-1),*transpose!=0,fgs_info[*handle],NULL,us);
    //

    // gs_gather_many_acc(u,u,vn,gsh->map_local[0^transpose],dom,op); 
    {
      for(k=0;k<vn;++k) {
	l_map = map;
	while((i=*l_map++)!=-(uint)1) {
	  t=ud[i+k*dstride];
	  j=*l_map++;
	  do { 
	    t += ud[j+k*dstride];
	  } while((j=*l_map++)!=-(uint)1);
	  ud[i+k*dstride]=t;
	}								
      }

    }
    // --
    if(dtrans==0) {
      // wtf?!?!
      // gs_init_many_acc(u,vn,gsh->flagged_primaries,dom,op);
      {
	for(k=0;k<vn;++k) {
	  l_map = map;
	  while((i=*l_map++)!=-(uint)1) {
	    ud[i+k*dstride]=0.0;
	  }
	}
      }
    }
    /* post receives */
    sbufp = pw_exec_recvs(bufp,vn*gs_dom_size[l_dom],comm,&pwd->comm[recv],pwd->req);
    /* fill send buffer */
    // gs_scatter_vec_to_many_acc(sendbuf,data,vn,pwd->map[send],dom);
    {
      q = sbufp;
      for(k=0;k<vn;k++) {
	l_map = map;
	while((i=*l_map++)!=-(uint)1) {				  
	  t=ud[i+k*dstride];
	  j=*l_map++;
	  do {
	    ((double*)q)[j*vn]=t;
	  } while((j=*l_map++)!=-(uint)1); 
	}								
	q+=gs_dom_size[l_dom];
      }
    }
    /* post sends */
    pw_exec_sends(sbufp,vn*gs_dom_size[l_dom],comm,&pwd->comm[send],&pwd->req[pwd->comm[recv].n]);
    comm_wait(pwd->req,pwd->comm[0].n+pwd->comm[1].n);
    /* gather using recv buffer */
    // gs_gather_vec_to_many_acc(data,buf,vn,pwd->map[recv],dom,op);
    {
      q = bufp;
      for(k=0;k<vn;k++) {
	l_map = map;
	while((i=*l_map++)!=-(uint)1) { 
	  t=ud[i+k*dstride];
	  j=*l_map++; 
	  do {
	    t += ((double*)q)[j*vn];
	  } while((j=*l_map++)!=-(uint)1);
	  ud[i+k*dstride]=t;
	}                                                                          
	q+=gs_dom_size[l_dom];
      }
    }
    // --
    // gs_scatter_many_acc(u,u,vn,gsh->map_local[1^transpose],dom); 
    {
      for(k=0;k<vn;++k) {
	l_map = map;
	while((i=*l_map++)!=-(uint)1) {
	  t=ud[i+k*dstride];
	  j=*l_map++; 
	  do {
	    ud[j+k*dstride]=t; 
	  } while((j=*l_map++)!=-(uint)1);
	}
      }
    }
  }
}
