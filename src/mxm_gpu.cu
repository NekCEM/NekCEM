/*
 * mxm_gpu.cu
 *  @author azamat, mmin
 *  @since  July 13, 2012
 */

#include <stdio.h>
#include <cuda.h>
#include <cuda_runtime.h>

#define KERNEL  1
#define TILE   16 //autotune-able
#define VERBOSE 0
#if VERBOSE
int dbg=1;
#else
int dbg=0;
#endif
static int once=0;
cudaEvent_t tstart, tstop, start, stop;
float kern=0.0f, xfer=0.0f;

#define onceMallocMemcpy(x,dbg) do{                                \
  if ((x)->sync&0x1) {                                             \
    cudaMalloc(&(x)->dev,(x)->sz);                                 \
    if(dbg){                                                       \
      printf("cudaMalloc'ed:     %s, %d B\n",(x)->vname,(x)->sz);  \
    }                                                              \
    (x)->sync^=0x1;                                                \
  }                                                                \
  if ((x)->sync&0x2) {                                             \
    cudaEventRecord(tstart,0);                                     \
    cudaMemcpy((x)->dev,(x)->host,(x)->sz,cudaMemcpyHostToDevice); \
    cudaEventRecord(tstop,0); cudaEventSynchronize(tstop);         \
    if(dbg){                                                       \
      cudaEventElapsedTime(&xfer,tstart,tstop);                    \
      printf("cudaMemcpy'ed H2D: %s, %d B, %f ms, %.2f MB/s\n",(x)->vname,(x)->sz,xfer,(1e3f*(x)->sz)/(xfer*(1<<20)));  \
    }                                                              \
    (x)->sync^=0x2;                                                \
  }                                                                \
}while(0)
#define onceMemcpyFree(x,dbg) do{                                  \
  if ((x)->sync&0x4) {                                             \
    cudaEventRecord(tstart,0);                                     \
    cudaMemcpy((x)->host,(x)->dev,(x)->sz,cudaMemcpyDeviceToHost); \
    cudaEventRecord(tstop,0); cudaEventSynchronize(tstop);         \
    if(dbg){                                                       \
      cudaEventElapsedTime(&xfer,tstart,tstop);                    \
      printf("cudaMemcpy'ed D2H: %s, %d B, %f ms, %.2f MB/s\n",(x)->vname,(x)->sz,xfer,(1e3f*(x)->sz)/(xfer*(1<<20)));  \
    }                                                              \
    (x)->sync^=0x4;                                                \
  }                                                                \
  if ((x)->sync&0x8) {                                             \
    cudaFree((x)->dev);                                            \
    if(dbg){                                                       \
      printf("cudaFree'ed:       %s\n",(x)->vname);                \
    }                                                              \
    (x)->sync^=0x8;                                                \
  }                                                                \
}while(0)

extern "C" {
  struct memptr {
    int sync; //sync flags: 0x1->allocate, 0x2->copy H2D, 0x4->copy D2H, 0x8->deallocate
    int sz;
    double* host;
    double* dev;
    char* vname;
  };
  typedef struct memptr memptr_t;
  void mxm_std_gpu_(double* a, int* m, double* b, int* n, double* c, int* p);
  void local_grad3_gpu_(
    memptr_t *u1r, memptr_t *u1s, memptr_t *u1t,
    memptr_t *u2r, memptr_t *u2s, memptr_t *u2t,
    memptr_t *u3r, memptr_t *u3s, memptr_t *u3t,
    memptr_t *u1 , memptr_t *u2 , memptr_t *u3 ,
    memptr_t *mp_d, memptr_t *mp_dt,
    int *n, int *nelts, int *lpts1, int *rank);
  void curl_gpu_(
    memptr_t *u1r, memptr_t *u1s, memptr_t *u1t,
    memptr_t *u2r, memptr_t *u2s, memptr_t *u2t,
    memptr_t *u3r, memptr_t *u3s, memptr_t *u3t,
    memptr_t *rxmn,memptr_t *sxmn,memptr_t *txmn,
    memptr_t *rymn,memptr_t *symn,memptr_t *tymn,
    memptr_t *rzmn,memptr_t *szmn,memptr_t *tzmn,
    memptr_t *w1,  memptr_t *w2,  memptr_t *w3,
    memptr_t *w3mn,
    int *nxyz, int *nelts, int *lpts1);
}


// basic curl kernel impl
__global__ void curl_vanilla(
    const double* __restrict__ rxmn,const double* __restrict__ rymn,const double* __restrict__ rzmn,
    const double* __restrict__ sxmn,const double* __restrict__ symn,const double* __restrict__ szmn,
    const double* __restrict__ txmn,const double* __restrict__ tymn,const double* __restrict__ tzmn,
    const double* __restrict__ u1r, const double* __restrict__ u1s, const double* __restrict__ u1t,
    const double* __restrict__ u2r, const double* __restrict__ u2s, const double* __restrict__ u2t,
    const double* __restrict__ u3r, const double* __restrict__ u3s, const double* __restrict__ u3t,
    const double* __restrict__ w3mn,const int nxyz, const int nelts,const int lpts1,
    double* __restrict__ w1
  ){
  const int tid=blockIdx.x*blockDim.x+threadIdx.x;
  double w3mk;
  int k=0;
  for(int e=0; e<nelts; e++){
    k=e*nxyz+tid;
    w3mk=w3mn[tid];

    w1[k]= w3mk*u3r[k]*rymn[k]
         + w3mk*u3s[k]*symn[k]
         + w3mk*u3t[k]*tymn[k]
         - w3mk*u2r[k]*rzmn[k]
         - w3mk*u2s[k]*szmn[k]
         - w3mk*u2t[k]*tzmn[k];

    w1[k+lpts1]
         = w3mk*u1r[k]*rzmn[k]
         + w3mk*u1s[k]*szmn[k]
         + w3mk*u1t[k]*tzmn[k]
         - w3mk*u3r[k]*rxmn[k]
         - w3mk*u3s[k]*sxmn[k]
         - w3mk*u3t[k]*txmn[k];

    w1[k+2*lpts1]
         = w3mk*u2r[k]*rxmn[k]
         + w3mk*u2s[k]*sxmn[k]
         + w3mk*u2t[k]*txmn[k]
         - w3mk*u1r[k]*rymn[k]
         - w3mk*u1s[k]*symn[k]
         - w3mk*u1t[k]*tymn[k];
  }
}

// basic multi-mxm impl
__global__ void mxm_vanilla(const double* __restrict__ a, const int m,
                            const double* __restrict__ b, const int n,
                            double* __restrict__ c, const int p,
                            const int nelts, const int ldims){
  const int row=blockIdx.y*blockDim.y+threadIdx.y;
  const int col=blockIdx.x*blockDim.x+threadIdx.x;
  if(row<m && col<p){ //eliminate out-of-bounds threads
    double s;
    int lda=( ldims&0x1)    *m*n    //if a's bit (0x1) is set, its leading dim is of size m*n 
      , ldb=((ldims&0x2)>>1)*n*p
      , ldc=((ldims&0x4)>>2)*m*p
      , ldi=((ldims&0x8)>>3)*m*n*p; //for inner dimensions
    if(ldims<8){ //no inner iterations
      for(int e=0; e<nelts; e++){
        s=0.0;
        for(int k=0; k<n; k++){
          s+=a[e*lda+k*m+row]*b[e*ldb+col*n+k];
        }
        c[e*ldc+col*m+row]=s;
      }
    }else{
      for(int e=0; e<nelts; e++){
        for(int i=0; i<m; i++){
          s=0.0;
          for(int k=0; k<n; k++){
            s+=a[e*ldi+i*lda+k*m+row]*b[col*n+k];
          }
          c[e*ldi+i*ldc+col*m+row]=s;
        }
      }
    }
  }
}


// mxm with 1D arrays
__global__ void mxm_1d(double* a, const int m, double* b, const int n, double* c, const int p){
  const int i=blockIdx.x*blockDim.x+threadIdx.x;
  if (i<m){
    for(int k=0; k<p; k++){
      double s=0.0;
      for(int j=0; j<n; j++){
        s+=a[j*m+i]*b[k*n+j];
      }
      c[k*m+i]=s;
    }
  }
}


// mxm with 2D arrays
__global__ void mxm_shared(double* a, const int m, double* b, const int n, double* c, const int p){
  __shared__ double as[TILE][TILE];
  __shared__ double bs[TILE][TILE];
  int bx=blockIdx.x, by=blockIdx.y, tx=threadIdx.x, ty=threadIdx.y;
  const int row=by*TILE+ty;
  const int col=bx*TILE+tx;
  double s=0.0;
  for(int t=0;t<m/TILE;t++){
    as[ty][tx]=a[col*m+t*TILE+tx];
    bs[ty][tx]=b[col*n+t*TILE+ty];
    __syncthreads();
    for(int k=0; k<TILE; k++){
      s+=as[ty][k]*bs[k][tx];
    }
    __syncthreads();
    c[col*m+row]=s;
  }
}


// globally-visible basic mxm implementation for small matrices
void mxm_std_gpu_(double* a, int* m, double* b, int* n, double* c, int* p){
  /*device variables*/
  double *dev_a, *dev_b, *dev_c;
  int sizeofA=*m*(*n)*sizeof(double)
    , sizeofB=*n*(*p)*sizeof(double)
    , sizeofC=*m*(*p)*sizeof(double);
  /*malloc and memcopy data H2D*/
  cudaMalloc(&dev_a,sizeofA);
  cudaMalloc(&dev_b,sizeofB);
  cudaMalloc(&dev_c,sizeofC);
  cudaMemcpy(dev_a,a,sizeofA,cudaMemcpyHostToDevice);
  cudaMemcpy(dev_b,b,sizeofB,cudaMemcpyHostToDevice);
  /*thread dimensions*/
  dim3 dimBlock, dimGrid;
#if KERNEL==1
  dimBlock.x=TILE; dimGrid.x=(*p+dimBlock.x-1)/dimBlock.x;
  dimBlock.y=TILE; dimGrid.y=(*m+dimBlock.y-1)/dimBlock.y;
  mxm_vanilla<<<dimGrid,dimBlock>>>(dev_a,*m,dev_b,*n,dev_c,*p,1,0);
#elif KERNEL==2
  dimBlock.x=TILE; dimGrid.x=(*m+dimBlock.x-1)/dimBlock.x;
  mxm_1d<<<dimGrid,dimBlock>>>(dev_a,*m,dev_b,*n,dev_c,*p);
#else
  dimBlock.x=TILE; dimGrid.x=(*p+dimBlock.x-1)/dimBlock.x;
  dimBlock.y=TILE; dimGrid.y=(*m+dimBlock.y-1)/dimBlock.y;
  mxm_shared<<<dimGrid,dimBlock>>>(dev_a,*m,dev_b,*n,dev_c,*p);
#endif
  /*memcopy D2H*/
  cudaMemcpy(c,dev_c,sizeofC,cudaMemcpyDeviceToHost);
  cudaFree(dev_a);
  cudaFree(dev_b);
  cudaFree(dev_c);
}


// sets up the aggregated mxm kernel launch
void mxm_gpu2(double* a, int as, int m
             ,double* b, int bs, int n
             ,double* c, int cs, int p
             ,int nelts, int mask, int dev){
  cudaSetDevice(dev);
  /*device variables*/
  double *dev_a, *dev_b, *dev_c;
  int sizeofA=as*sizeof(double)
    , sizeofB=bs*sizeof(double)
    , sizeofC=cs*sizeof(double);
  /*malloc and memcopy H2D*/
  cudaMalloc(&dev_a,sizeofA);
  cudaMalloc(&dev_b,sizeofB);
  cudaMalloc(&dev_c,sizeofC);
  cudaMemcpy(dev_a,a,sizeofA,cudaMemcpyHostToDevice);
  cudaMemcpy(dev_b,b,sizeofB,cudaMemcpyHostToDevice);
  /*thread grid dimensions*/
  dim3 dimBlock, dimGrid;
  dimBlock.x=TILE; dimGrid.x=(p+dimBlock.x-1)/dimBlock.x;
  dimBlock.y=TILE; dimGrid.y=(m+dimBlock.y-1)/dimBlock.y;
  mxm_vanilla<<<dimGrid,dimBlock>>>(dev_a,m, dev_b,n, dev_c,p, nelts,mask);
  /*memcopy D2H*/
  cudaMemcpy(c,dev_c,sizeofC,cudaMemcpyDeviceToHost);
  cudaFree(dev_a);
  cudaFree(dev_b);
  cudaFree(dev_c);
}

//=============================================================================
// sets up the aggregated mxm kernel launch
void mxm_gpu_agg(memptr_t *a, int m
                ,memptr_t *b, int n
                ,memptr_t *c, int p
                ,int nelts, int mask, int dev){
  cudaSetDevice(dev);
  /*malloc and memcopy H2D*/
  onceMallocMemcpy(a,dbg);
  onceMallocMemcpy(b,dbg);
  onceMallocMemcpy(c,dbg);
  /*thread grid dimensions*/
  dim3 dimBlock, dimGrid;
  dimBlock.x=TILE; dimGrid.x=(p+dimBlock.x-1)/dimBlock.x;
  dimBlock.y=TILE; dimGrid.y=(m+dimBlock.y-1)/dimBlock.y;
  mxm_vanilla<<<dimGrid,dimBlock>>>(a->dev,m, b->dev,n, c->dev,p, nelts,mask);
  /*memcopy D2H and dealloc*/
  onceMemcpyFree(a,dbg);
  onceMemcpyFree(b,dbg);
  onceMemcpyFree(c,dbg);
}


//=============================================================================
/**
 * Performs aggregated mxm for all elements at once.
 *
 * foreach e in 0..nelts
 *   u@r_{NxN^2} = d_{NxN} * u@_{NxN^2}^{e} // here @ is either 1, 2 or 3
 *   foreach k in 0..N
 *     u@s_{NxN}^{k} = u@_{NxN}^{k,e} * dt_{NxN}
 *   u@t_{N^2xN} = u@_{N^2xN}^{e} * dt_{NxN}
 */
void local_grad3_gpu_(memptr_t *u1r, memptr_t *u1s, memptr_t *u1t,  
                      memptr_t *u2r, memptr_t *u2s, memptr_t *u2t,  
                      memptr_t *u3r, memptr_t *u3s, memptr_t *u3t,  
                      memptr_t *u1 , memptr_t *u2 , memptr_t *u3 ,  
                      memptr_t *d,   memptr_t *dt,
                      int *n, int *nelts, int *lpts1, int *rank){
  int n2=*n*(*n);
  if (!once) {
    d->vname   = "d";
    dt->vname  = "dt";
    u1r->vname = "u1r";
    u1s->vname = "u1s";
    u1t->vname = "u1t";
    u2r->vname = "u2r";
    u2s->vname = "u2s";
    u2t->vname = "u2t";
    u3r->vname = "u3r";
    u3s->vname = "u3s";
    u3t->vname = "u3t";
    u1->vname  = "u1";
    u2->vname  = "u2";
    u3->vname  = "u3";
    cudaEventCreate(&tstart); cudaEventCreate(&tstop);
    cudaEventCreate(&start);  cudaEventCreate(&stop);
  }

  // select the device
  int devs = 0;
  cudaGetDeviceCount(&devs);
  int devid;
  if (devs==1) {
    devid = 0;
  } else {
    devid = *rank%2;
  }
  cudaSetDevice(devid);

  onceMallocMemcpy(d,  dbg);
  // u1,u2,u3 are contiguous, do a single transfer
  u1->sz=*lpts1*3*sizeof(double);
  onceMallocMemcpy(u1, dbg);
  u2->dev=u1->dev+(*lpts1);
  u3->dev=u2->dev+(*lpts1);
  onceMallocMemcpy(u1r,dbg);
  onceMallocMemcpy(u2r,dbg);
  onceMallocMemcpy(u3r,dbg);

  /*thread grid dimensions*/
  dim3 dimBlock, dimGrid;
  dimBlock.x=TILE; dimBlock.y=TILE;

  //         d_{NxN}   * u*_{NxN^2}= u*r_{NxN^2}  foreach e
  dimGrid.x=(n2+dimBlock.x-1)/dimBlock.x;
  dimGrid.y=(*n+dimBlock.y-1)/dimBlock.y;
  cudaEventRecord(start,0);
  mxm_vanilla<<<dimGrid,dimBlock>>>(d->dev,*n,  u1->dev,*n, u1r->dev,n2, *nelts,6);
  mxm_vanilla<<<dimGrid,dimBlock>>>(d->dev,*n,  u2->dev,*n, u2r->dev,n2, *nelts,6);
  mxm_vanilla<<<dimGrid,dimBlock>>>(d->dev,*n,  u3->dev,*n, u3r->dev,n2, *nelts,6);
  cudaEventRecord(stop,0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&kern,start,stop);
  if(dbg){
    printf("r kernel time:     %f ms\n",kern);
  }

  //         u*_{NxN}  * dt_{NxN}  = u*s_{NxN}    foreach e,k
  onceMallocMemcpy(dt, dbg);
  onceMallocMemcpy(u1s,dbg);
  onceMallocMemcpy(u2s,dbg);
  onceMallocMemcpy(u3s,dbg);
  dimGrid.x=(*n+dimBlock.x-1)/dimBlock.x;
  dimGrid.y=(*n+dimBlock.y-1)/dimBlock.y;
  cudaEventRecord(start,0);
  mxm_vanilla<<<dimGrid,dimBlock>>>(u1->dev,*n, dt->dev,*n, u1s->dev,*n, *nelts,13);
  mxm_vanilla<<<dimGrid,dimBlock>>>(u2->dev,*n, dt->dev,*n, u2s->dev,*n, *nelts,13);
  mxm_vanilla<<<dimGrid,dimBlock>>>(u3->dev,*n, dt->dev,*n, u3s->dev,*n, *nelts,13);
  cudaEventRecord(stop,0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&kern,start,stop);
  if(dbg){
    printf("s kernel time:     %f ms\n",kern);
  }

  //         u*_{N^2xN}* dt_{NxN}  = u*t_{N^2xN}  foreach e
  onceMallocMemcpy(u1t,dbg);
  onceMallocMemcpy(u2t,dbg);
  onceMallocMemcpy(u3t,dbg);
  dimGrid.x=(*n+dimBlock.x-1)/dimBlock.x;
  dimGrid.y=(n2+dimBlock.y-1)/dimBlock.y;
  cudaEventRecord(start,0);
  mxm_vanilla<<<dimGrid,dimBlock>>>(u1->dev,n2, dt->dev,*n, u1t->dev,*n, *nelts,5);
  mxm_vanilla<<<dimGrid,dimBlock>>>(u2->dev,n2, dt->dev,*n, u2t->dev,*n, *nelts,5);
  mxm_vanilla<<<dimGrid,dimBlock>>>(u3->dev,n2, dt->dev,*n, u3t->dev,*n, *nelts,5);
  cudaEventRecord(stop,0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&kern,start,stop);
  if(dbg){
    printf("t kernel time:     %f ms\n",kern);
  }
  // nothing to copy D2H or to free
  cudaDeviceSynchronize();
}

//=============================================================================
// Sets up the curl kernel
void curl_gpu_(memptr_t *u1r,  memptr_t *u1s,  memptr_t *u1t,
               memptr_t *u2r,  memptr_t *u2s,  memptr_t *u2t,
               memptr_t *u3r,  memptr_t *u3s,  memptr_t *u3t,
               memptr_t *rxmn, memptr_t *sxmn, memptr_t *txmn,
               memptr_t *rymn, memptr_t *symn, memptr_t *tymn,
               memptr_t *rzmn, memptr_t *szmn, memptr_t *tzmn,
               memptr_t *w1,   memptr_t *w2,   memptr_t *w3,
               memptr_t *w3mn, int *nxyz, int *nelts, int *lpts1){
  if (!once){
    rxmn->vname="rxmn"; sxmn->vname="sxmn"; txmn->vname="txmn";
    rymn->vname="rymn"; symn->vname="symn"; tymn->vname="tymn";
    rzmn->vname="rzmn"; szmn->vname="szmn"; tzmn->vname="tzmn";
    w3mn->vname="w3mn";
    w1->vname="w1"; w2->vname="w2"; w3->vname="w3";
    once=1;
  }
  /*malloc and memcopy H2D*/
  onceMallocMemcpy(rxmn,dbg);
  onceMallocMemcpy(rymn,dbg);
  onceMallocMemcpy(rzmn,dbg);
  onceMallocMemcpy(sxmn,dbg);
  onceMallocMemcpy(symn,dbg);
  onceMallocMemcpy(szmn,dbg);
  onceMallocMemcpy(txmn,dbg);
  onceMallocMemcpy(tymn,dbg);
  onceMallocMemcpy(tzmn,dbg);
  onceMallocMemcpy(w3mn,dbg);
  onceMallocMemcpy(u1r, dbg);
  onceMallocMemcpy(u1s, dbg);
  onceMallocMemcpy(u1t, dbg);
  onceMallocMemcpy(u2r, dbg);
  onceMallocMemcpy(u2s, dbg);
  onceMallocMemcpy(u2t, dbg);
  onceMallocMemcpy(u3r, dbg);
  onceMallocMemcpy(u3s, dbg);
  onceMallocMemcpy(u3t, dbg);
  // w1,w2,w3 are contiguous, do a single transfer
  w1->sz=*lpts1*3*sizeof(double);
  onceMallocMemcpy(w1,  dbg);
  /*thread grid dimensions*/
  dim3 dimBlock, dimGrid;
  dimBlock.x=*nxyz; dimGrid.x=(15+dimBlock.x-1)/dimBlock.x;
  curl_vanilla<<<dimGrid,dimBlock>>>(
    rxmn->dev,rymn->dev,rzmn->dev,
    sxmn->dev,symn->dev,szmn->dev,
    txmn->dev,tymn->dev,tzmn->dev,
    u1r->dev, u1s->dev, u1t->dev,
    u2r->dev, u2s->dev, u2t->dev,
    u3r->dev, u3s->dev, u3t->dev,
    w3mn->dev,*nxyz,*nelts, *lpts1,
    w1->dev
  );
  onceMemcpyFree(rxmn,dbg);
  onceMemcpyFree(rymn,dbg);
  onceMemcpyFree(rzmn,dbg);
  onceMemcpyFree(sxmn,dbg);
  onceMemcpyFree(symn,dbg);
  onceMemcpyFree(szmn,dbg);
  onceMemcpyFree(txmn,dbg);
  onceMemcpyFree(tymn,dbg);
  onceMemcpyFree(tzmn,dbg);
  onceMemcpyFree(w3mn,dbg);
  onceMemcpyFree(u1r, dbg);
  onceMemcpyFree(u1s, dbg);
  onceMemcpyFree(u1t, dbg);
  onceMemcpyFree(u2r, dbg);
  onceMemcpyFree(u2s, dbg);
  onceMemcpyFree(u2t, dbg);
  onceMemcpyFree(u3r, dbg);
  onceMemcpyFree(u3s, dbg);
  onceMemcpyFree(u3t, dbg);
  onceMemcpyFree(w1,  dbg);
}

