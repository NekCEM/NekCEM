#include <openacc.h>

void fgs_fields_acc(const sint *handle,
		    double *u, const sint *stride, const sint *n,
		    const sint *dom, const sint *op, const sint *transpose)
{
  struct pw_data *pwd;
  struct comm    *comm;
  buffer         *buf;
  const unsigned recv = 0^*transpose, send = 1^*transpose;
  uint    i,j,k,bs,uds,dstride,dtrans,vn,m_size;
  int    *l_map,*map,*t_map,*fp_map,*snd_map,*rcv_map;
  char   *bufp,*sbufp,*q;
  double  *dbufp,*dsbufp,*dq;
  double  t,*ud;
  gs_dom  l_dom;

  //  bs = vn*gs_dom_size[l_dom]*(fgs_info[*handle]->r.buffer_size);
  buf = &static_buffer;
  bs = (*n)*8*(fgs_info[*handle]->r.buffer_size);
  buffer_reserve(buf,bs);
  fgs_check_parms(*handle,*dom,*op,"gs_op_fields",__LINE__);
  if(*n<0) return;

  uds = (*stride) * (*n); // Get size of u in number of doubles
  pwd     = fgs_info[*handle]->r.data;
  comm    = &fgs_info[*handle]->comm;
  m_size  = fgs_info[*handle]->handle_size;
  map     = (int*)(fgs_info[*handle]->map_local[0^*transpose]);
  t_map   = (int*)(fgs_info[*handle]->map_local[1^*transpose]);
  fp_map  = (int*)(fgs_info[*handle]->flagged_primaries);
  snd_map = (int*)(pwd->map[send]);
  rcv_map = (int*)(pwd->map[recv]);

#pragma enter data present_or_create(ud[uds],dbufp[bs],dsbufp[bs],l_map[m_size]) 
  {
#pragma enter data present_or_copyin(t_map[m_size],map[m_size],fp_map[m_size],snd_map[m_size],rcv_map[m_size],buf[bs])
    {
#pragma data present(u[uds])
      {
  //  acc_present(u,uds);

  //  acc_pcreate(ud,uds);
  //  acc_pcreate(dbufp,bs);
  //  acc_pcreate(dsbufp,bs);
  //  acc_pcreate(buf,bs);
  //  acc_pcopyin(t_map,m_size);
  //  acc_pcopyin(map,m_size);
  //  acc_pcopyin(fp_map,m_size);//(int*)(fgs_info[*handle]->flagged_primaries);
  //  acc_pcopyin(snd_map,m_size);//(int*)(pwd->map[send]);
  //  acc_pcopyin(rcv_map,m_size);//(int*)(pwd->map[recv]);
  //  acc_pcreate(l_map,m_size);


    // Get the map and other pointers we need.
#pragma acc kernels present(u[uds],ud[uds])
  {
  // Setup a "double" version of u
    ud = u;                 // Get pointer from u
    dstride = *stride;
    //    l_dom   = fgs_dom[*dom];

    dtrans  = *transpose;
    vn      = *n;
    
    // Create temp buffer for gather/scatter and send/recv
    dbufp = (double*)buf->ptr;
    //    dbufp = (double*)bufp;
  }


  //#pragma acc data create(bufp[0:bs]) present(ud[0:uds])
  {
    // The below implementing cgs_many()/gs_aux():
    //
    // gs_aux_acc(u,mode_many,dn,fgs_dom[*dom],(gs_op_t)(*op-1),*transpose!=0,fgs_info[*handle],NULL,us);
    //

    // gs_gather_many_acc(u,u,vn,gsh->map_local[0^transpose],dom,op); 
    {
#pragma acc parallel loop gang vector present(ud[0:uds],l_map[0:m_size],map[0:m_size])
      for(k=0;k<vn;++k) {
	l_map = map;
	while((i=*l_map++)!=-(uint)1) {
	  t=ud[i+k*dstride];
	  j=*l_map++;
	  do { 
	    //            printf("gather i: %d j: %d t: %lf out: %lf \n",i,j,t,ud[j+k*dstride]);
	    t += ud[j+k*dstride];
	    //            printf("gather2 i: %d j: %d t: %lf out: %lf \n",i,j,t,ud[j+k*dstride]);
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
#pragma acc parallel loop gang vector present(ud[0:uds],l_map[0:m_size],fp_map[0:m_size])
	for(k=0;k<vn;++k) {
          l_map = fp_map;
	  while((i=*l_map++)!=-(uint)1) {
	    ud[i+k*dstride]=0.0;
	  }
	}
      }
    }

    /* post receives */
    dsbufp = (double*)pw_exec_recvs((char*)dbufp,vn*8,comm,&pwd->comm[recv],pwd->req);

    //    printf("bufp: %lX sbu %lX\n",bufp,sbufp);
    /* fill send buffer */
    // gs_scatter_many_to_vec_acc(sendbuf,data,vn,pwd->map[send],dom);
    {
      //      dq = dsbufp;
#pragma acc parallel loop gang vector present(ud[0:uds],l_map[0:m_size],snd_map[0:m_size],dsbufp[0:bs])
      for(k=0;k<vn;k++) {
	l_map = snd_map;
	while((i=*l_map++)!=-(uint)1) {				  
	  t=ud[i+k*dstride];
	  j=*l_map++;
	  do {
	    ///            printf("scatter_##T i: %d j: %d t: %lf out: %lf : %d\n",i,j,t,((double*)q)[j*vn]);
	    dsbufp[j*vn]=t;
	    //            printf("scatter_##T2 i: %d j: %d t: %lf out: %lf : %d\n",i,j,t,((double*)q)[j*vn]);
	  } while((j=*l_map++)!=-(uint)1); 
	}		
	//        printf("p1: %lX\n",q);
	dsbufp++;//sizeof(double)=gs_dom_size[l_dom];
	//        printf("p1: %lX\n",q);
      }
    }
    /* post sends */

#pragma update host(dsbufp)
    pw_exec_sends((char*)dsbufp,vn*8,comm,&pwd->comm[send],&pwd->req[pwd->comm[recv].n]);
    comm_wait(pwd->req,pwd->comm[0].n+pwd->comm[1].n);
#pragma update device(dbufp) 
    /* gather using recv buffer */
    // gs_gather_vec_to_many_acc(data,buf,vn,pwd->map[recv],dom,op);
    {
      //      dq = dbufp;
#pragma acc parallel loop gang vector present(ud[0:uds],l_map[0:m_size],rcv_map[0:m_size],dbufp[bs])
      for(k=0;k<vn;k++) {
	l_map = rcv_map;
	while((i=*l_map++)!=-(uint)1) { 
	  t=ud[i+k*dstride];
	  j=*l_map++; 
	  do {
	    //            printf("gather i: %d j: %d t: %lf out: %lf \n",i,j,t,((double*)q)[j*vn]);
	    t += dbufp[j*vn];
	    //            printf("gather2 i: %d j: %d t: %lf out: %lf \n",i,j,t,((double*)q)[j*vn]);
	  } while((j=*l_map++)!=-(uint)1);
	  ud[i+k*dstride]=t;
	  //          printf("out: %lf\n",ud[i+k*dstride]);
	}                               
	//        printf("p2: %lX\n",q);
	dbufp++;//sizeof(double)=gs_dom_size[l_dom];
	//        printf("p2: %lX\n",q);
      }
    }
    // --
    // gs_scatter_many_acc(u,u,vn,gsh->map_local[1^transpose],dom); 
    {
#pragma acc parallel loop gang vector present(ud[0:uds],l_map[0:m_size],t_map[0:m_size])
      for(k=0;k<vn;++k) {
	l_map = t_map;
	while((i=*l_map++)!=-(uint)1) {
	  t=ud[i+k*dstride];
	  j=*l_map++; 
	  do {
	    //            printf("scatter_##T i: %d j: %d t: %lf out: %lf index: %d\n",i,j,t,ud[j+k*dstride],j+k*dstride);
	    ud[j+k*dstride]=t; 
	    //            printf("scatter_##T2 i: %d j: %d t: %lf out: %lf\n",i,j,t,ud[j+k*dstride]);
	  } while((j=*l_map++)!=-(uint)1);
	}
      }
    }
  }
      }}}
}
