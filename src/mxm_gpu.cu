#include <stdio.h>
#include <cuda.h>
#include <cuda_runtime.h>

#define KERNEL 1
#define TILE 16

extern "C" {
  void local_grad3_gpu_(double* u1r, double* u1s, double* u1t,  
                        double* u2r, double* u2s, double* u2t,  
                        double* u3r, double* u3s, double* u3t,  
                        double* u1 , double* u2 , double* u3 ,  
                        double* dxm1,   int* n,      int* nelts);
  void mxm_std_gpu_(double* a, int* m, double* b, int* n, double* c, int* p);
}


void print_array(double* a, int m, int n){
  int i,j,k=0;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      printf("array[%d][%d]=%E\n",i,j,a[k++]);
    }
  }
}
__global__ void mxm_vanilla(double* a, const int m, double* b, const int n, double* c, const int p
                           ,const int nelts, const int ldims){
  const int row=blockIdx.y*blockDim.y+threadIdx.y;
  const int col=blockIdx.x*blockDim.x+threadIdx.x;
  if(row<m && col<p){//eliminate out-of-bounds threads
    int lda=(ldims&0x1)*m*n //if a's bit (0x1) is set, its leading dim is of size m*n 
      , ldb=(ldims&0x2)*n*p
      , ldc=(ldims&0x4)*m*p;
    for(int e=0; e<nelts; e++){ // might need to launch 1 thread per element
      double s=0.0;
      for(int k=0; k<n; k++){
        s+=a[e*lda+k*m+row]*b[e*ldb+col*n+k];
      }
      c[e*ldc+col*m+row]=s;
    }
  }
}
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
void mxm_std_gpu_(double* a, int* m, double* b, int* n, double* c, int* p){
  //printf("mxm_gpu: m=%d,n=%d,p=%d\n",*m,*n,*p);
  //print_array(c,*m,*p);
  /*device variables*/
  double *dev_a, *dev_b, *dev_c;
  int sizeofA=*m*(*n)*sizeof(double)
    , sizeofB=*n*(*p)*sizeof(double)
    , sizeofC=*m*(*p)*sizeof(double);
  /*malloc and memcopy data from host to device*/
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
  //printf("mxm_gpu: dimGrid.x=%d,dimGrid.y=%d\n",dimGrid.x,dimGrid.y);
  /*memcopy from device to host*/
  cudaMemcpy(c,dev_c,sizeofC,cudaMemcpyDeviceToHost);
  cudaFree(dev_a);
  cudaFree(dev_b);
  cudaFree(dev_c);
  cudaDeviceSynchronize();
}
void mxm_gpu2(double* a, int as, int m
             ,double* b, int bs, int n
             ,double* c, int cs, int p
             ,int nelts, int mask){
  //printf("mxm_gpu: m=%d,n=%d,p=%d\n",*m,*n,*p);
  //print_array(c,*m,*p);
  /*device variables*/
  double *dev_a, *dev_b, *dev_c;
  int sizeofA=as*sizeof(double)
    , sizeofB=bs*sizeof(double)
    , sizeofC=cs*sizeof(double);
  /*malloc and memcopy data from host to device*/
  cudaMalloc(&dev_a,sizeofA);
  cudaMalloc(&dev_b,sizeofB);
  cudaMalloc(&dev_c,sizeofC);
  cudaMemcpy(dev_a,a,sizeofA,cudaMemcpyHostToDevice);
  cudaMemcpy(dev_b,b,sizeofB,cudaMemcpyHostToDevice);
  /*thread dimensions*/
  dim3 dimBlock, dimGrid;
#if KERNEL==1
  dimBlock.x=TILE; dimGrid.x=(p+dimBlock.x-1)/dimBlock.x;
  dimBlock.y=TILE; dimGrid.y=(m+dimBlock.y-1)/dimBlock.y;
  mxm_vanilla<<<dimGrid,dimBlock>>>(dev_a,m, dev_b,n, dev_c,p, nelts,mask);
#elif KERNEL==2
  dimBlock.x=TILE; dimGrid.x=(m+dimBlock.x-1)/dimBlock.x;
  mxm_1d<<<dimGrid,dimBlock>>>(dev_a,m,dev_b,n,dev_c,p);
#else
  dimBlock.x=TILE; dimGrid.x=(p+dimBlock.x-1)/dimBlock.x;
  dimBlock.y=TILE; dimGrid.y=(m+dimBlock.y-1)/dimBlock.y;
  mxm_shared<<<dimGrid,dimBlock>>>(dev_a,m,dev_b,n,dev_c,p);
#endif
  //printf("mxm_gpu: dimGrid.x=%d,dimGrid.y=%d\n",dimGrid.x,dimGrid.y);
  /*memcopy from device to host*/
  cudaMemcpy(c,dev_c,sizeofC,cudaMemcpyDeviceToHost);
  cudaFree(dev_a);
  cudaFree(dev_b);
  cudaFree(dev_c);
  cudaDeviceSynchronize();
}
void local_grad3_gpu_(double* u1r, double* u1s, double* u1t,  
                      double* u2r, double* u2s, double* u2t,  
                      double* u3r, double* u3s, double* u3t,  
                      double* u1 , double* u2 , double* u3 ,  
                      double* dxm1,   int* n  ,    int* nelts){
  // foreach e in 0..nelts
  //   u*r_{NxN^2} = d_{NxN} * u*_{NxN^2}^{e} // * is either 1, 2 or 3
  //   foreach k in 0..N
  //     u*s_{NxN}^{k} = u*_{NxN}^{k,e} * dt_{NxN}
  //   u*t_{N^2xN} = u*_{N^2xN}^{e} * dt_{NxN}
//  int n2=*n*(*n)
//    , n3=*n*n2;
//    , npts=n3*(*nelts);
//  // calc u1r
//  mxm_gpu2(dxm1,n2  ,*n
//          ,u1  ,npts,*n
//          ,u1r ,npts,n2
//          ,*nelts,2);
}

