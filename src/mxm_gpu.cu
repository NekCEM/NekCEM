#include <stdio.h>
#include <cuda.h>
#include <cuda_runtime.h>

#define KERNEL 2
#define TILE 16

extern "C" {
  void local_grad3_gpu_(double* u1r, double* u1s, double* u1t,  
                        double* u2r, double* u2s, double* u2t,  
                        double* u3r, double* u3s, double* u3t,  
                        double* u1 , double* u2 , double* u3 ,  
       int* nxyz, int* nelt, int* npts, double* dxm1, int* N);
}

extern "C" {
  void mxm_gpu_(double* a, int* m, double* b, int* n, double* c, int* p);
}
void print_array(double* a, int m, int n){
  int i,j,k=0;
  for(j=0; j<n; j++){
    for(i=0; i<m; i++){
      printf("array[%d][%d]=%E\n",i,j,a[k++]);
    }
  }
}
__global__ void mxm_vanilla(double* a, const int m, double* b, const int n, double* c, const int p){
  const int row=blockIdx.y*blockDim.y+threadIdx.y;
  const int col=blockIdx.x*blockDim.x+threadIdx.x;
  double s=0.0;
  if (row<m && col<p){
    for(int k=0; k<n; k++){
      s+=a[row+k*m]*b[k+col*n];
    }
    c[row+col*m]=s;
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
      //printf("%d.%d:s=%E\n",blockIdx.x,threadIdx.x,s);
    }
  }
}
void mxm_gpu_(double* a, int* m, double* b, int* n, double* c, int* p){
  //printf("mxm_gpu: m=%d,n=%d,p=%d\n",*m,*n,*p);
  //print_array(c,*m,*p);
  /*device variables*/
  double *dev_a, *dev_b, *dev_c;
  int sizeofA=*m*(*n)*sizeof(double);
  int sizeofB=*n*(*p)*sizeof(double);
  int sizeofC=*m*(*p)*sizeof(double);
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
  mxm_vanilla<<<dimGrid,dimBlock>>>(dev_a,*m,dev_b,*n,dev_c,*p);
#else
  dimBlock.x=TILE; dimGrid.x=(*m+dimBlock.x-1)/dimBlock.x;
  mxm_1d<<<dimGrid,dimBlock>>>(dev_a,*m,dev_b,*n,dev_c,*p);
#endif
  //printf("mxm_gpu: dimGrid.x=%d,dimGrid.y=%d\n",dimGrid.x,dimGrid.y);
  /*memcopy from device to host*/
  cudaMemcpy(c,dev_c,sizeofC,cudaMemcpyDeviceToHost);
  cudaFree(dev_a);
  cudaFree(dev_b);
  cudaFree(dev_c);
  cudaDeviceSynchronize();
}
