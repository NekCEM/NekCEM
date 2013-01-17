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

extern "C" {
  void local_grad3_gpu_(double* u1r, double* u1s, double* u1t,  
                        double* u2r, double* u2s, double* u2t,  
                        double* u3r, double* u3s, double* u3t,  
                        double* u1 , double* u2 , double* u3 ,  
                        double* dxm, double* dxtm, int* n, int* nelts, int* rank);
  void mxm_std_gpu_(double* a, int* m, double* b, int* n, double* c, int* p);
}


// basic multi-mxm impl
__global__ void mxm_vanilla(double* a, const int m, double* b, const int n, double* c, const int p
                           ,const int nelts, const int ldims){
  const int row=blockIdx.y*blockDim.y+threadIdx.y;
  const int col=blockIdx.x*blockDim.x+threadIdx.x;
  if(row<m && col<p){ //eliminate out-of-bounds threads
    if(nelts==1){ //single mxm
      double s=0.0;
      for(int k=0; k<n; k++){
        s+=a[k*m+row]*b[col*n+k];
      }
      c[col*m+row]=s;
    }else{ //multiple mxm's
      int lda=(ldims&0x1)*m*n //if a's bit (0x1) is set, its leading dim is of size m*n 
        , ldb=(ldims&0x2)*n*p
        , ldc=(ldims&0x4)*m*p
        , ldai=(ldims&0x8)*m*n //for a's inner dimension
        , ldci=(ldims&0x8)*m*p;
      for(int e=0; e<nelts; e++){ // might need to launch 1 thread per element
        if(ldims<8){ //no inner iterations
          double s=0.0;
          for(int k=0; k<n; k++){
            s+=a[e*lda+k*m+row]*b[e*ldb+col*n+k];
          }
          c[e*ldc+col*m+row]=s;
        }else{
          for(int i=0; i<m; i++){
            double s=0.0;
            for(int k=0; k<n; k++){
              s+=a[e*lda+i*ldai+k*m+row]*b[e*ldb+col*n+k];
            }
            c[e*ldc+i*ldci+col*m+row]=s;
          }
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
                      double* d  , double* dt , int* n, int* nelts, int* rank){
  int n2=*n*(*n), n3=*n*n2, npts=n3*(*nelts);

  int devs = 0;
  cudaGetDeviceCount(&devs);
  int devid = *rank%2;
  if (devs==1) {
    //       d_{NxN}   *  u*_{NxN^2} = u*r_{NxN^2}   foreach e
    mxm_gpu2(d,n2,*n,     u1,npts,*n,  u1r,npts,n2,  *nelts,6, 0);
    mxm_gpu2(d,n2,*n,     u2,npts,*n,  u2r,npts,n2,  *nelts,6, 0);
    mxm_gpu2(d,n2,*n,     u3,npts,*n,  u3r,npts,n2,  *nelts,6, 0);
  
    //       u*_{NxN}  *  dt_{NxN}  =  u*s_{NxN}     foreach e,k
    mxm_gpu2(u1,npts,*n,  dt,n2,*n,    u1s,npts,*n,  *nelts,13, 0);
    mxm_gpu2(u2,npts,*n,  dt,n2,*n,    u2s,npts,*n,  *nelts,13, 0);
    mxm_gpu2(u3,npts,*n,  dt,n2,*n,    u3s,npts,*n,  *nelts,13, 0);
  
    //       u*_{N^2xN} * dt_{NxN}  =  u*t_{N^2xN}   foreach e
    mxm_gpu2(u1,npts,n2,  dt,n2,*n,    u1t,npts,*n,  *nelts,5, 0);
    mxm_gpu2(u2,npts,n2,  dt,n2,*n,    u2t,npts,*n,  *nelts,5, 0);
    mxm_gpu2(u3,npts,n2,  dt,n2,*n,    u3t,npts,*n,  *nelts,5, 0);
  } else {
    // todo: fork threads or do async launches
    //       d_{NxN}   *  u*_{NxN^2} = u*r_{NxN^2}   foreach e
    mxm_gpu2(d,n2,*n,     u1,npts,*n,  u1r,npts,n2,  *nelts,6, devid);
    mxm_gpu2(d,n2,*n,     u2,npts,*n,  u2r,npts,n2,  *nelts,6, devid);
    mxm_gpu2(d,n2,*n,     u3,npts,*n,  u3r,npts,n2,  *nelts,6, devid);
  
    //       u*_{NxN}  *  dt_{NxN}  =  u*s_{NxN}     foreach e,k
    mxm_gpu2(u1,npts,*n,  dt,n2,*n,    u1s,npts,*n,  *nelts,13, devid);
    mxm_gpu2(u2,npts,*n,  dt,n2,*n,    u2s,npts,*n,  *nelts,13, devid);
    mxm_gpu2(u3,npts,*n,  dt,n2,*n,    u3s,npts,*n,  *nelts,13, devid);
  
    //       u*_{N^2xN} * dt_{NxN}  =  u*t_{N^2xN}   foreach e
    mxm_gpu2(u1,npts,n2,  dt,n2,*n,    u1t,npts,*n,  *nelts,5, devid);
    mxm_gpu2(u2,npts,n2,  dt,n2,*n,    u2t,npts,*n,  *nelts,5, devid);
    mxm_gpu2(u3,npts,n2,  dt,n2,*n,    u3t,npts,*n,  *nelts,5, devid);
  }
}

