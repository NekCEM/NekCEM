#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "c99.h"
#include "name.h"
#include "fail.h"
#include "types.h"
#include "comm.h"
#include "mem.h"
#include "gs_defs.h"
#include "gs.h"

typedef double T;
const gs_dom dom = gs_double;

static void test(const struct comm *comm)
{
  struct gs_data *gsh;
  const uint np = comm->np;
  slong *id = tmalloc(slong,np+4);
  T *v = tmalloc(T,np+4);
  uint i;
  id[0] = -(slong)(np+10+3*comm->id);
  for(i=0;i<np;++i) id[i+1] = -(sint)(i+1);
  id[np+1] = comm->id+1;
  id[np+2] = comm->id+1;
  id[np+3] = np-comm->id;
  gsh = gs_setup(id,np+4,comm,0,gs_auto,1);
  free(id);
  
  for(i=0;i<np+4;++i) v[i] = 1;
  gs(v,dom,gs_add,0,gsh,0);
  if(comm->id==0) for(i=0;i<np+4;++i) printf("%g\n",v[i]);
  if(comm->id==0) printf("\n");
  for(i=0;i<np+4;++i) v[i] = 1;
  gs(v,dom,gs_add,1,gsh,0);
  if(comm->id==0) for(i=0;i<np+4;++i) printf("%g\n",v[i]);

  gs_free(gsh);
  free(v);
}

int main(int narg, char *arg[])
{
  comm_ext world; int np;
  FILE *inp;
  struct comm comm;
  int *targetCoreIndexing,index;
  int *sendbuf,*recvbuf,*sendcounts,*displs;
  int ret,totalElements,maxNp,v1,v2,v3,v4;
  int arrayOffset;
  int maxElementsPerCore,elementIndex,targetCore,i,j,nid,k;
  char buffer[1024];
  char inpFile[128];
  int **mat;


  printf("Beginning\n");
#ifdef MPI
  MPI_Init(&narg,&arg);
  world = MPI_COMM_WORLD;
  MPI_Comm_size(world,&np);
#else
  world=0, np=1;
#endif
  

  printf("After mpi\n");
  comm_init(&comm,world);
  MPI_Comm_rank(world,&nid);



  if(nid==0){
    inp = fopen(arg[1],"r");
    mat = (int **) malloc(sizeof(int *)*np);
    //Only supports 2d right now
    fgets(buffer,1024,inp);
    sscanf(buffer,"%d%d%d%d%d%d%d",&totalElements,&v1,&v2,&maxNp,&i,&j,&k);
    //scanf("%d%d%d%d%d%d%d",&totalElements,&v1,&v2,&maxNp,&i,&j,&k);
    //    fscanf(inp,"%d%d%d%d%d%d%d",&totalElements,&v1,&v2,&maxNp,&i,&j,&k);
    maxElementsPerCore = ceil((double)totalElements);
    printf("maxElementsPerCore: %d\n",maxElementsPerCore);
    MPI_Bcast(&maxElementsPerCore,1,MPI_INT,0,world);
    //Allocate memory for each core
    //Allocate memory for the elements
    for(i=0;i<np;i++) mat[i] = (int*) malloc(sizeof(int)*4*maxElementsPerCore);

    //Set an array for indexing the elements in the matrix
    targetCoreIndexing = (int *) malloc(sizeof(int)*np);
    for(i=0;i<np;i++) {
      targetCoreIndexing[i] = 0;
    }

    //read in the data
    for(i=0;i<totalElements;i++){

      fgets(buffer,1024,inp);
      sscanf(buffer,"%d%d%d%d%d", &elementIndex, &v1, &v2, &v3, &v4);
      //scanf("%d%d%d%d%d", &elementIndex, &v1, &v2, &v3, &v4);
      //      ret = fscanf(inp,"%d%d%d%d%d", &elementIndex, &v1, &v2, &v3, &v4);
      //if(ret == EOF) printf("EOF\n");

      targetCore = elementIndex%np;

      arrayOffset = targetCoreIndexing[targetCore];

      //Store the vertices

      k = 0+4*arrayOffset;
      mat[targetCore][k] = v1;

      k++;
      mat[targetCore][k] = v2;

      k++;
      mat[targetCore][k] = v3;

      k++;
      mat[targetCore][k] = v4;

      targetCoreIndexing[targetCore] = targetCoreIndexing[targetCore]+1;
    }
    
    sendbuf    = malloc(sizeof(int)*4*maxElementsPerCore*np);
    sendcounts = malloc(sizeof(int)*np);
    displs     = malloc(sizeof(int)*np);
    k=0;
    for(i=0;i<np;i++) {
      sendcounts[i] = 4*maxElementsPerCore;
      displs[i]     = i*4*maxElementsPerCore;
      //Fill the send buffer
      for(j=0;j<targetCoreIndexing[i];j++){
        sendbuf[k] = mat[i][0+4*j];
        printf("sendbuf[%d]: %d\n",k,sendbuf[k]);
        k++;
        sendbuf[k] = mat[i][1+4*j];
        printf("sendbuf[%d]: %d\n",k,sendbuf[k]);
        k++;
        sendbuf[k] = mat[i][2+4*j];
        printf("sendbuf[%d]: %d\n",k,sendbuf[k]);
        k++;
        sendbuf[k] = mat[i][3+4*j];
        printf("sendbuf[%d]: %d\n",k,sendbuf[k]);
      }
    }
    recvbuf = malloc(sizeof(int)*4*maxElementsPerCore);
    MPI_Scatterv(sendbuf,sendcounts,displs,MPI_INT,recvbuf,4*maxElementsPerCore,MPI_INT,0,world);
  } else {
    MPI_Bcast(&maxElementsPerCore,1,MPI_INT,0,world);
    printf("Nid: %d maxEle: %d\n",nid,maxElementsPerCore);
    recvbuf = malloc(sizeof(int)*4*maxElementsPerCore);
    MPI_Scatterv(sendbuf,sendcounts,displs,MPI_INT,recvbuf,4*maxElementsPerCore,MPI_INT,0,world);
  }

  for(i=0;i<maxElementsPerCore;i++){
    printf("%d: ",nid);
    k = i;
    printf("%d ",recvbuf[k]);
    k++;
    printf("%d ",recvbuf[k]);
    k++;
    printf("%d ",recvbuf[k]);
    k++;
    printf("%d \n",recvbuf[k]);
  }

  //  test(&comm);
  
  comm_free(&comm);
  if(nid==0){
    for(i=0;i<np;i++) free(mat[i]);
    free(mat);
    free(targetCoreIndexing);
    free(sendbuf);
    free(sendcounts);
    free(displs);
  }
  free(recvbuf);

#ifdef MPI
  MPI_Barrier(world);
  MPI_Finalize();
#endif

  return 0;
}
