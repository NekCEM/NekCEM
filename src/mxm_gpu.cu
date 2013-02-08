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

extern "C" {
  struct memptr {
    int sync;
    int sz;
    double* host;
    double* dev;
    char* vname;
  };
  typedef struct memptr memptr_t;
  void mxm_std_gpu_(double* a, int* m, double* b, int* n, double* c, int* p);
  void local_grad3_gpu_(
    double* u1r, double* u1s, double* u1t,
    double* u2r, double* u2s, double* u2t,
    double* u3r, double* u3s, double* u3t,
    double* u1 , double* u2 , double* u3 ,
    memptr_t *mp_d, memptr_t *mp_dt, int* n, int* nelts, int* rank);
  void curl_gpu_(
    double* u1r, double* u1s, double* u1t,
    double* u2r, double* u2s, double* u2t,
    double* u3r, double* u3s, double* u3t,
    double* rxmn,double* sxmn,double* txmn,
    double* rymn,double* symn,double* tymn,
    double* rzmn,double* szmn,double* tzmn,
    double* w1,  double* w2,  double* w3, double* w3m, int* nxyz, int* nelts);
}


// basic curl kernel impl
__global__ void curl_vanilla(
    double* rxmn,double* rymn,double* rzmn,
    double* sxmn,double* symn,double* szmn,
    double* txmn,double* tymn,double* tzmn,
    double* u1r, double* u1s, double* u1t,
    double* u2r, double* u2s, double* u2t,
    double* u3r, double* u3s, double* u3t,
    double* w3m, const int nxyz, const int nelts,
    double* w1,  double* w2,  double* w3){
  const int tid=blockIdx.x*blockDim.x+threadIdx.x;
  double w3mk;
  int k=0;
  for(int e=0; e<nelts; e++){
    k=e*nxyz+tid;
    w3mk=w3m[k];

    w1[k]= w3mk*u3r[k]*rymn[k]
         + w3mk*u3s[k]*symn[k]
         + w3mk*u3t[k]*tymn[k]
         - w3mk*u2r[k]*rzmn[k]
         - w3mk*u2s[k]*szmn[k]
         - w3mk*u2t[k]*tzmn[k];

    w2[k]= w3mk*u1r[k]*rzmn[k]
         + w3mk*u1s[k]*szmn[k]
         + w3mk*u1t[k]*tzmn[k]
         - w3mk*u3r[k]*rxmn[k]
         - w3mk*u3s[k]*sxmn[k]
         - w3mk*u3t[k]*txmn[k];

    w3[k]= w3mk*u2r[k]*rxmn[k]
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
  cudaDeviceSynchronize();
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
  cudaDeviceSynchronize();
}

//=============================================================================
// sets up the aggregated mxm kernel launch
void mxm_gpu_agg(memptr_t *a, int m
                ,memptr_t *b, int n
                ,memptr_t *c, int p
                ,int nelts, int mask, int dev){
  cudaSetDevice(dev);
  /*malloc and memcopy H2D*/
  if (a->sync&0x1) { // if need to malloc
    cudaMalloc(&a->dev,a->sz);
    #if(VERBOSE)
    printf("cudaMalloc'ed: %s\n",a->vname);
    #endif
    a->sync^=0x1;
  }
  if (b->sync&0x1) {
    cudaMalloc(&b->dev,b->sz);
    #if(VERBOSE)
    printf("cudaMalloc'ed: %s\n",b->vname);
    #endif
    b->sync^=0x1;
  }
  if (c->sync&0x1) {
    cudaMalloc(&c->dev,c->sz);
    #if(VERBOSE)
    printf("cudaMalloc'ed: %s\n",c->vname);
    #endif
    c->sync^=0x1;
  }
  if (a->sync&0x2) { // if need to memcpy H2D
    cudaMemcpy(a->dev,a->host,a->sz,cudaMemcpyHostToDevice);
    #if(VERBOSE)
    printf("cudaMemcpy'ed H2D: %s\n",a->vname);
    #endif
    a->sync^=0x2;
  }
  if (b->sync&0x2) {
    cudaMemcpy(b->dev,b->host,b->sz,cudaMemcpyHostToDevice);
    #if(VERBOSE)
    printf("cudaMemcpy'ed H2D: %s\n",b->vname);
    #endif
    b->sync^=0x2;
  }
  if (c->sync&0x2) {
    cudaMemcpy(c->dev,c->host,c->sz,cudaMemcpyHostToDevice);
    #if(VERBOSE)
    printf("cudaMemcpy'ed H2D: %s\n",c->vname);
    #endif
    c->sync^=0x2;
  }
  /*thread grid dimensions*/
  dim3 dimBlock, dimGrid;
  dimBlock.x=TILE; dimGrid.x=(p+dimBlock.x-1)/dimBlock.x;
  dimBlock.y=TILE; dimGrid.y=(m+dimBlock.y-1)/dimBlock.y;
  mxm_vanilla<<<dimGrid,dimBlock>>>(a->dev,m, b->dev,n, c->dev,p, nelts,mask);
  /*memcopy D2H and dealloc*/
  if (a->sync&0x4) { // if need to memcpy D2H
    cudaMemcpy(a->host,a->dev,a->sz,cudaMemcpyDeviceToHost);
    #if(VERBOSE)
    printf("cudaMemcpy'ed D2H: %s\n",a->vname);
    #endif
    a->sync^=0x4;
  }
  if (b->sync&0x4) {
    cudaMemcpy(b->host,b->dev,b->sz,cudaMemcpyDeviceToHost);
    #if(VERBOSE)
    printf("cudaMemcpy'ed D2H: %s\n",b->vname);
    #endif
    b->sync^=0x4;
  }
  if (c->sync&0x4) {
    cudaMemcpy(c->host,c->dev,c->sz,cudaMemcpyDeviceToHost);
    #if(VERBOSE)
    printf("cudaMemcpy'ed D2H: %s\n",c->vname);
    #endif
    c->sync^=0x4;
  }
  if (a->sync&0x8) { // if need to dealloc
    cudaFree(a->dev);
    #if(VERBOSE)
    printf("cudaFree'ed: %s\n",a->vname);
    #endif
    a->sync^=0x8;
  }
  if (b->sync&0x8) {
    cudaFree(b->dev);
    #if(VERBOSE)
    printf("cudaFree'ed: %s\n",b->vname);
    #endif
    b->sync^=0x8;
  }
  if (c->sync&0x8) {
    cudaFree(c->dev);
    #if(VERBOSE)
    printf("cudaFree'ed: %s\n",c->vname);
    #endif
    c->sync^=0x8;
  }
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
void local_grad3_gpu_(double* u1r, double* u1s, double* u1t,  
                      double* u2r, double* u2s, double* u2t,  
                      double* u3r, double* u3s, double* u3t,  
                      double* u1 , double* u2 , double* u3 ,  
                      memptr_t *mp_d, memptr_t *mp_dt, int* n, int* nelts, int* rank){
  int n2=*n*(*n), n3=*n*n2, npts=n3*(*nelts);
  int szdbl = sizeof(double);
  // device pointers
  double *du1 =NULL, *du2 =NULL, *du3 =NULL;
  double *du1r=NULL, *du2r=NULL, *du3r=NULL;
  double *du1s=NULL, *du2s=NULL, *du3s=NULL;
  double *du1t=NULL, *du2t=NULL, *du3t=NULL;
  // sync flags: 0x1->allocate, 0x2->copy H2D, 0x4->copy D2H, 0x8->deallocate
  mp_d->vname = "d";
  mp_dt->vname = "dt";
  memptr_t mp_u1  = {3, szdbl*npts, u1, du1 , "u1" };
  memptr_t mp_u2  = {3, szdbl*npts, u2, du2 , "u2" };
  memptr_t mp_u3  = {3, szdbl*npts, u3, du3 , "u3" };
  memptr_t mp_u1r = {13, szdbl*npts, u1r,du1r, "u1r"};
  memptr_t mp_u2r = {13, szdbl*npts, u2r,du2r, "u2r"};
  memptr_t mp_u3r = {13, szdbl*npts, u3r,du3r, "u3r"};
  memptr_t mp_u1s = {13, szdbl*npts, u1s,du1s, "u1s"};
  memptr_t mp_u2s = {13, szdbl*npts, u2s,du2s, "u2s"};
  memptr_t mp_u3s = {13, szdbl*npts, u3s,du3s, "u3s"};
  memptr_t mp_u1t = {13, szdbl*npts, u1t,du1t, "u1t"};
  memptr_t mp_u2t = {13, szdbl*npts, u2t,du2t, "u2t"};
  memptr_t mp_u3t = {13, szdbl*npts, u3t,du3t, "u3t"};

  // select the device
  int devs = 0;
  cudaGetDeviceCount(&devs);
  int devid;
  if (devs==1) {
    devid = 0;
  } else {
    devid = *rank%2;
  }

  // todo: fork threads or do async launches
  //         d_{NxN}   * u*_{NxN^2}= u*r_{NxN^2}  foreach e
  mxm_gpu_agg(mp_d,*n,   &mp_u1,*n,  &mp_u1r,n2,  *nelts,6, devid);
  mxm_gpu_agg(mp_d,*n,   &mp_u2,*n,  &mp_u2r,n2,  *nelts,6, devid);
  mxm_gpu_agg(mp_d,*n,   &mp_u3,*n,  &mp_u3r,n2,  *nelts,6, devid);
  
  //         u*_{NxN}  * dt_{NxN}  = u*s_{NxN}    foreach e,k
  mxm_gpu_agg(&mp_u1,*n, mp_dt,*n,   &mp_u1s,*n,  *nelts,13, devid);
  mxm_gpu_agg(&mp_u2,*n, mp_dt,*n,   &mp_u2s,*n,  *nelts,13, devid);
  mxm_gpu_agg(&mp_u3,*n, mp_dt,*n,   &mp_u3s,*n,  *nelts,13, devid);
 
  //         u*_{N^2xN}* dt_{NxN}  = u*t_{N^2xN}  foreach e
  mp_u1.sync=8;
  mxm_gpu_agg(&mp_u1,n2, mp_dt,*n,   &mp_u1t,*n,  *nelts,5, devid);
  mp_u2.sync=8;
  mxm_gpu_agg(&mp_u2,n2, mp_dt,*n,   &mp_u2t,*n,  *nelts,5, devid);
  mp_u3.sync=8;
  mxm_gpu_agg(&mp_u3,n2, mp_dt,*n,   &mp_u3t,*n,  *nelts,5, devid);
}

//=============================================================================
// Sets up the curl kernel
void curl_gpu_(double* u1r, double* u1s, double* u1t,
               double* u2r, double* u2s, double* u2t,
               double* u3r, double* u3s, double* u3t,
               double* rxmn,double* sxmn,double* txmn,
               double* rymn,double* symn,double* tymn,
               double* rzmn,double* szmn,double* tzmn,
               double* w1,  double* w2,  double* w3, double* w3m, int* nxyz, int* nelts){
  /*device variables*/
  double *dev_rxmn, *dev_rymn, *dev_rzmn
        ,*dev_sxmn, *dev_symn, *dev_szmn
        ,*dev_txmn, *dev_tymn, *dev_tzmn
        ,*dev_u1r, *dev_u1s, *dev_u1t
        ,*dev_u2r, *dev_u2s, *dev_u2t
        ,*dev_u3r, *dev_u3s, *dev_u3t
        ,*dev_w1, *dev_w2, *dev_w3, *dev_w3m;
  int nptsz=*nxyz*(*nelts)*sizeof(double);
  /*malloc and memcopy H2D*/
  cudaMalloc(&dev_rxmn,nptsz);
  cudaMalloc(&dev_rymn,nptsz);
  cudaMalloc(&dev_rzmn,nptsz);
  cudaMalloc(&dev_sxmn,nptsz);
  cudaMalloc(&dev_symn,nptsz);
  cudaMalloc(&dev_szmn,nptsz);
  cudaMalloc(&dev_txmn,nptsz);
  cudaMalloc(&dev_tymn,nptsz);
  cudaMalloc(&dev_tzmn,nptsz);
  cudaMalloc(&dev_u1r, nptsz);
  cudaMalloc(&dev_u1s, nptsz);
  cudaMalloc(&dev_u1t, nptsz);
  cudaMalloc(&dev_u2r, nptsz);
  cudaMalloc(&dev_u2s, nptsz);
  cudaMalloc(&dev_u2t, nptsz);
  cudaMalloc(&dev_u3r, nptsz);
  cudaMalloc(&dev_u3s, nptsz);
  cudaMalloc(&dev_u3t, nptsz);
  cudaMalloc(&dev_w3m, nptsz);
  cudaMalloc(&dev_w1,  nptsz);
  cudaMalloc(&dev_w2,  nptsz);
  cudaMalloc(&dev_w3,  nptsz);
  cudaMemcpy(dev_rxmn,rxmn,nptsz,cudaMemcpyHostToDevice);
  cudaMemcpy(dev_rymn,rymn,nptsz,cudaMemcpyHostToDevice);
  cudaMemcpy(dev_rzmn,rzmn,nptsz,cudaMemcpyHostToDevice);
  cudaMemcpy(dev_sxmn,sxmn,nptsz,cudaMemcpyHostToDevice);
  cudaMemcpy(dev_symn,symn,nptsz,cudaMemcpyHostToDevice);
  cudaMemcpy(dev_szmn,szmn,nptsz,cudaMemcpyHostToDevice);
  cudaMemcpy(dev_txmn,txmn,nptsz,cudaMemcpyHostToDevice);
  cudaMemcpy(dev_tymn,tymn,nptsz,cudaMemcpyHostToDevice);
  cudaMemcpy(dev_tzmn,tzmn,nptsz,cudaMemcpyHostToDevice);
  cudaMemcpy(dev_u1r, u1r, nptsz,cudaMemcpyHostToDevice);
  cudaMemcpy(dev_u1s, u1s, nptsz,cudaMemcpyHostToDevice);
  cudaMemcpy(dev_u1t, u1t, nptsz,cudaMemcpyHostToDevice);
  cudaMemcpy(dev_u2r, u2r, nptsz,cudaMemcpyHostToDevice);
  cudaMemcpy(dev_u2s, u2s, nptsz,cudaMemcpyHostToDevice);
  cudaMemcpy(dev_u2t, u2t, nptsz,cudaMemcpyHostToDevice);
  cudaMemcpy(dev_u3r, u3r, nptsz,cudaMemcpyHostToDevice);
  cudaMemcpy(dev_u3s, u3s, nptsz,cudaMemcpyHostToDevice);
  cudaMemcpy(dev_u3t, u3t, nptsz,cudaMemcpyHostToDevice);
  cudaMemcpy(dev_w3m, w3m, nptsz,cudaMemcpyHostToDevice);
  /*thread grid dimensions*/
  dim3 dimBlock, dimGrid;
  dimBlock.x=*nxyz; dimGrid.x=(15+dimBlock.x-1)/dimBlock.x;
  curl_vanilla<<<dimGrid,dimBlock>>>(
    dev_rxmn,dev_rymn,dev_rzmn,
    dev_sxmn,dev_symn,dev_szmn,
    dev_txmn,dev_tymn,dev_tzmn,
    dev_u1r,dev_u1s,dev_u1t,
    dev_u2r,dev_u2s,dev_u2t,
    dev_u3r,dev_u3s,dev_u3t,
    dev_w3m,*nxyz,*nelts,
    dev_w1, dev_w2, dev_w3
  );
  cudaMemcpy(w1,dev_w1,nptsz,cudaMemcpyDeviceToHost);
  cudaMemcpy(w2,dev_w2,nptsz,cudaMemcpyDeviceToHost);
  cudaMemcpy(w3,dev_w3,nptsz,cudaMemcpyDeviceToHost);
  cudaFree(dev_rxmn);
  cudaFree(dev_rymn);
  cudaFree(dev_rzmn);
  cudaFree(dev_sxmn);
  cudaFree(dev_symn);
  cudaFree(dev_szmn);
  cudaFree(dev_txmn);
  cudaFree(dev_tymn);
  cudaFree(dev_tzmn);
  cudaFree(dev_u1r);
  cudaFree(dev_u1s);
  cudaFree(dev_u1t);
  cudaFree(dev_u2r);
  cudaFree(dev_u2s);
  cudaFree(dev_u2t);
  cudaFree(dev_u3r);
  cudaFree(dev_u3s);
  cudaFree(dev_u3t);
  cudaFree(dev_w3m);
  cudaFree(dev_w1);
  cudaFree(dev_w2);
  cudaFree(dev_w3);
  cudaDeviceSynchronize();
}

