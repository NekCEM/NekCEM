/********************************sparse_matrix.c*******************************
Module Name: sparse_matrix.c
Module Info:

Author:  Henry M. Tufo III

e-mail:  hmt@@cs.brown.edu (lgrind ... so excuse the double <at>)

sn-mail: Division of Applied Mathematics, 
	 Brown University,
         Box F
	 Providence, RI 02912

Tel:	 (401) 863-7666


Last Modification: 12.21.02
*********************************sparse_matrix.c******************************/


/********************************sparse_matrix.c*******************************
NOTES ON USAGE: 

*********************************sparse_matrix.c******************************/


/********************************sparse_matrix.c*******************************
FILE FORMAT: 
------------------------------ Begin File -------------------------------------

------------------------------ End   File -------------------------------------

Note: 
*********************************sparse_matrix.c******************************/

/* C modules for I/O etc. */
#if   defined NXSRC
#ifndef DELTA
#include <nx.h>
#endif

#elif defined MPISRC
#include <mpi.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include <unistd.h>

/* mine : const before types! */
#include "const.h"
#include "types.h"
#include "error.h"
#include "ivec.h"
#include "bss_malloc.h"
#include "adj_list.h"
#include "queue.h"
#include "stack.h"

/* struct for lanczos */
struct lanczos_node {
  double d;
  double e;
  double *v;
};

/* sparse matrix template */
typedef struct id{
  int lda, sda, nnz;
  char ty[80];
  char desc[80];
  int  *permutation;
  int  gs_handle;
  adj_node_ptr adj_list;
  int     *row_lens;
  int    **row_indices;
  double **row_vals;
  int     *col_lens;
  int    **col_indices;
  double **col_vals;
} matrix_id;
typedef struct id  *matrix_id_ptr;

/* struct for recursive bisection */
struct tree_node{
  int depth;
  int rank;
  int num_act;
  int *active;
  int num_sep;
  int *sep;
  double *rhs;
  struct tree_node *lc, *rc;
  int lb, ub;
};

struct rsb_info{
  int max_depth;
  int depth;
  struct tree_node *root;
  int nel;
  int nv;
  int nvu;
  int nvo;
  int n;
  int *vertex;
  int *vertex_map;  
  queue_ADT done_q;
  queue_ADT half_done_q;
  queue_ADT bad_q;
};


/* defines */
#define   MAX_LEN       101             /* max input line length        */
#define   MAX_LDA       100000          /* max size of matrix           */
#define   FLD_WTH       6               /* printing info                */
#define   MAX_MAT       100
#define   HEADER_LEN    4
#define   SPARSE        0.25
#define   TRI           1
#define   BI            2
#define CG_LIMIT       10000
#define CG_MR_LIMIT    250
#define INIT_CG_LIMIT  10
#define INIT_RES_LIMIT 1.0e-02 /*pow(2.0,-20.0)*/
#define RES_LIMIT      1.0e-14 /*pow(2.0,-20.0)*/
#define RES_ORTHOG     1.0e-05 /*pow(2.0,-20.0)*/
#define PI             acos(-1.0)

/* for dvec sorting */
#define   DSORT_OPT	6     
#define   DSORT_STACK	50000


static void *doffset_stack[2*DSORT_STACK];
static int   dsize_stack[DSORT_STACK];
static PTRINT pdsize_stack[DSORT_STACK];


#ifdef NOT
/* c defined function proto. Definition in respective .h file */
/* some conflicts amongst include files ... so comment out   */
/* Note: return values for fprintf and printf are ignored ... */
double sqrt(double x);
double fabs(double x);
FILE *fopen(char *name,char *mode);
int fclose(FILE *fp);
int fprintf(FILE * stream,const char *format, ...);
int printf(const char *format, ...);
char *strtok(char *s,const char *ct);
char *fgets(char * s,int n,FILE *stream);
void rewind(FILE * stream);
void *malloc(size_t size);
int atoi(const char *s);
void free(void *p);
void  exit(int status);
#endif

/* sparse matrix prototypes - most should be static? */
int log_2_floor(int n);
void compress_adj_list(void);
matrix_id_ptr extract_sub_matrix(int *m,int n);
struct tree_node *new_tree_node();
int SMI_csc(int *ptrs, int *cols, double *vals, int lda, int type);
int SMI_csc_map(int *map, int *ptrs, int *cols, double *vals, int lda, 
		int type);
double *lanczos(double *rhs, int n, int *ia, int *ja, double *vals);
double lap_sub_matrix(int *list, int n, int **ia, int **ja, double **vals);
void hmt_ortho(double *r, int n);
double *Init_guess(int n);
double *Init_rhs(int n);
double hmt_ddot(double *x1, double *x2, int n);
void rev_daxpy(double *x1, double scale, double *x2, int n);
void hmt_daxpy(double *x1, double scale, double *x2, int n);
void matrix_mult(double *res, double *v, int n, int *ia, int *ja, double *aij);
void skew_matrix_mult(double *Ap_k, double *p_k, int n, int *ia, int *ja, 
		      double *vals, double lambda);
double *cg_min_res(double *rhs, int n, int *ia, int *ja, double *vals, 
		   double lambda);
void tri_calc_evec(double *d, double *sd, int n, double lambda, double *ret);
int trid_slv(double *d, double *sd, int n, double lambda, double *b);
double *trid_mult(double *d, double *sd, int n, double *b);
void normalize(double *r, int n);
double median(double *fiedler, int n);
double *hmt_perturb(double *rhs, int n);
void check_row_sum(int n, int *ia, int *ja, double *vals);
int *separate_graph(int *list, int n, double *rhs);
void dump_laplacian(int n, int *ia, int *ja, double *vals);

/* dvec prototypes */
static void dvec_dump(double *v, int n, int tag, int tag2, char * s);
static void dvec_zero(double *arg1, int n);
static void dvec_copy(double *arg1, double *arg2, int n);
static void dvec_set(double *arg1, double arg2, int n);
static void dvec_sort(double *ar, int Size);
static void dvec_sort_companion(double *ar, int *ar2, int Size);

/* Global Variables (not visible) */
static matrix_id_ptr matrix_list[MAX_MAT];
static int num_mat=0;

/* hack - should make obsolete */
static matrix_id_ptr active_matrix;
static double start_nnz;
static double wt_subg_level;
static double wt_cuts_level;
static double sep_sz_level;
static int nnz_level;
static double start_nnz;
static int num_cuts, sda, *pairs;
static double *vals;
static int max_depth;
static int *color_proc, *color_sep;
static queue_ADT done_q=NULL;
static queue_ADT bad_q=NULL;
static queue_ADT half_done_q=NULL;
struct tree_node *root=NULL;
struct rsb_info *rsb=NULL;

/* included from rsb_driver.h */
#define _2D                  0
#define _3D                  1
static int max_repeats=1;
static double res_limit=1e-14;
static int rsb_dim=2;
static double *xc=NULL, *yc=NULL, *zc=NULL;
static int *pre_nek_color;

#ifdef NOT
/* mistake to try jumping into fortran fp routines */
double fddot(int n, double *x, int incx, double *y, int incy);
void fdaxpy(long int n, double da, double *dx, long int incx, 
	   double *dy, long int incy);


double 
fddot(int n, double *x, int incx, double *y, int incy)
{
  ddot (&n, x, &incx, y, &incy);
}

void 
fdaxpy(long int n, double da, double *dx, long int incx, 
	   double *dy, long int incy)
{
  daxpy (&n, &da, dx, &incx, dy, &incy);
}
#endif


/********************************sparse_matrix.c*******************************
Function: 

Input : 
Output: 
Return: 
Description:  
*********************************sparse_matrix.c******************************/
double *
Init_rhs_random(double *rhs, int n)
{
  int i;
  double fact, *dptr1;


  /* slight perturbation? */  
  dvec_zero(rhs,n);
  fact = 2.0/sqrt(n);
  dptr1 = rhs;

  srand((int) time(NULL));
  for (i=0; i<n; i++)
    {*dptr1++ = fact*(((RAND_MAX-rand())/((double)RAND_MAX))-0.5);}
  return(rhs);
}


/********************************sparse_matrix.c*******************************
Function: 

Input : 
Output: 
Return: 
Description:  
*********************************sparse_matrix.c******************************/
double *
hmt_perturb(double *rhs, int n)
{
  int i;
  double fact, *dptr1;

  /* slight perturbation? */  
  normalize(rhs,n);

  dptr1 = rhs;
  for (i=0; i<n; i++)
    {
      fact = *dptr1/100.0;
      *dptr1++ += fact*(((RAND_MAX-rand())/((double)RAND_MAX))-0.5);
    }
  return(rhs);
}


/********************************sparse_matrix.c*******************************
Function: Init_rhs_coord()

Input : 
Output: 
Return: 
Description:  

This routine finds the x,y,z diameters and divides about the smallest. Note 
that for this to work properly the original mesh from nekton must be "square"
w.r.t. the coordiate system. Note also assumes that bw_...() isn't used on
the initial system!

Warning: relying in geometric info can be hazardous to your health 
and habit forming!!!
*********************************sparse_matrix.c******************************/
double *Init_rhs_coord(double *rhs, int *list, int n)
{
  int i, who;
  double x_max, x_min, y_max, y_min, z_max, z_min;
  double dx, dy, dz, dmin;
  double scale;
  double *tmp;


  if (rsb_dim==_2D)
    {
      if (!xc||!yc)
	{error_msg_fatal("Hey ... no coordinate info avail for 2D Init_rhs_coord!\n");}
    }
  else
    {
      if (!xc||!yc||!zc)
	{error_msg_fatal("Hey ... no coordinate info avail for 3D Init_rhs_coord!\n");}
    }

  tmp=rhs;

  x_max=y_max=z_max=DBL_MIN;
  x_min=y_min=z_min=DBL_MAX;

  if (rsb_dim==_2D)
    {
      for (i=0; i<n; i++)
	{
	  x_max=MAX(x_max,xc[i]);
	  x_min=MIN(x_min,xc[i]);
	  y_max=MAX(y_max,yc[i]);
	  y_min=MIN(y_min,yc[i]);
	}

      dx = x_max-x_min;
      dy = y_max-y_min;

      i=0;
      if (dx==0.0) {dx=FLT_MIN; i++;}
      if (dy==0.0) {dy=FLT_MIN; i++;}

      if ((i==2)&&(n>1))
	{error_msg_fatal("A bit anorexic are we?\n");}

      dmin=MAX(dx,dy);

      if (dx==dmin)
	{who=X;}
      else if (dy==dmin)
	{who=Y;}
      else
	{error_msg_fatal("Init_rhs_coord() :: 2D no min?\n");}

      scale = 1.0/sqrt((double) n);

      switch (who) {
      case X:
	dmin = dx/2.0 + x_min;
	for (i=0; i<n; i++)
	  {
	    if (xc[i]<dmin)
	      {*tmp++ = -scale;}
	    else
	      {*tmp++ =  scale;}
	  }
	break;
      case Y:
	dmin = dy/2.0 + y_min;
	for (i=0; i<n; i++)
	  {
	    if (yc[i]<dmin)
	      {*tmp++ = -scale;}
	    else
	      {*tmp++ =  scale;}
	  }
	break;
      default:
	error_msg_fatal("Init_rhs_coord() :: 2D bad axis?\n");
	break;
      }
    }
  else
    {
      for (i=0; i<n; i++)
	{
	  x_max=MAX(x_max,xc[i]);
	  x_min=MIN(x_min,xc[i]);
	  y_max=MAX(y_max,yc[i]);
	  y_min=MIN(y_min,yc[i]);
	  z_max=MAX(z_max,zc[i]);
	  z_min=MIN(z_min,zc[i]);
	}
      dx = x_max-x_min;
      dy = y_max-y_min;
      dz = z_max-z_min;

      i=0;
      if (dx==0.0) {dx=FLT_MIN; i++;}
      if (dy==0.0) {dy=FLT_MIN; i++;}
      if (dz==0.0) {dz=FLT_MIN; i++;}

      if ((i==3)&&(n>1))
	{error_msg_fatal("A bit anorexic are we?\n");}

      dmin=MAX(dx,dy);
      dmin=MAX(dmin,dz);

      if (dx==dmin)
	{who=X;}
      else if (dy==dmin)
	{who=Y;}
      else if (dz==dmin)
	{who=Z;}
      else
	{error_msg_fatal("Init_rhs_coord() :: 3D no min?\n");}

      scale = 1.0/sqrt((double) n);

#ifdef DEBUG
      printf("scale = %g\n",scale);
      printf("who=%d, dmin=%g\n",who,dmin);
      printf("dx=%g, dy=%g, dz=%g\n",dx,dy,dz);
#endif

      switch (who) {
      case X:
	for (i=0; i<n; i++)
	  {
	    dmin = dx/2.0 + x_min;
	    if (xc[i]<dmin)
	      {*tmp++ = -scale;}
	    else
	      {*tmp++ =  scale;}
	  }
	break;
      case Y:
	for (i=0; i<n; i++)
	  {
	    dmin = dy/2.0 + y_min;
	    if (yc[i]<dmin)
	      {*tmp++ = -scale;}
	    else
	      {*tmp++ =  scale;}
	  }
	break;
      case Z:
	for (i=0; i<n; i++)
	  {
	    dmin = dz/2.0 + z_min;
	    if (zc[i]<dmin)
	      {*tmp++ = -scale;}
	    else
	      {*tmp++ =  scale;}
	  }
	break;
      default:
	error_msg_fatal("Init_rhs_coord() :: 2D bad axis?\n");
	break;
      }
    }

#ifdef DEBUG
  error_msg_warning("Bad RHS!!!\n");
  error_msg_warning("RHS to follow ... \n");
  for (i=0; i<n; i++)
    {printf("%g\n",rhs[i]);}
#endif

  return(rhs);
}



/********************************sparse_matrix.c*******************************
Function: lap_sub_matrix()

Input : 
Output: 
Return: 
Description:  
*********************************sparse_matrix.c******************************/
double
lap_sub_matrix(int *list, register int n, int **ia, int **ja, double **vals)
{
  register int i, j, k, elm;
  register int *iptr_ia, *iptr_ja;
  register int *iptr_row, *iptr_list;
  double *dptr_vals, *v_hold, row_sum;
  int **rows;
  double **v;
  double *dptr, wt;
  int ct;

#ifdef DEBUG
  error_msg_warning("lap_sub_matrix() :: begin\n");
#endif

  /* need len and indice information */
  rows = active_matrix->row_indices;
  v    = active_matrix->row_vals;

  /* get space for ia offsets and initialize to all zeros */
  *ia = iptr_ia = (int *) bss_malloc(INT_LEN*(n+2));
  ivec_zero(iptr_ia,n+2);

  /* get space for ja and set ia */
  *iptr_ia = 0;
  for (i=j=0; i<n; i++)
    {
      iptr_list = list;
      iptr_row = *(rows + *(iptr_list+i));

      k=0;
      ct=0;
      while ((elm=*iptr_row++)>=0)
	{
	  while ((*iptr_list>=0)&&(elm > *iptr_list))
	    {iptr_list++;}
	  if (*iptr_list<0) break;
	  if (elm == *iptr_list)
	    {k++; iptr_list++;}
	}

#ifdef DEBUG
      if ((k<2) && (n>1))
	{error_msg_warning("lap_sub_matrix() :: disconnected graph?\n");}
#endif

      j+=k;
      *(iptr_ia + 1) = *(iptr_ia) + k;
      iptr_ia++;
    }
  *(iptr_ia + 1) = -1;

  /* ok ... get space for ja and vals */
  *ja = iptr_ja = (int *) bss_malloc(INT_LEN*j);
  ivec_zero(iptr_ja,j);
  *vals = dptr = (double *) bss_malloc(sizeof(double)*j);

  /* finish laplacian map ==> (ia,ja,vals) */
  iptr_ia = *ia;
  dvec_zero(dptr,j);
  
  for (wt=0.0,i=0; i<n; i++)
    {
      iptr_list = list;
      k = *(iptr_list+i);
      iptr_row = *(rows + k);
      dptr_vals = *(v + k);

      j=0;
      row_sum = 0.0;
      while ((elm=*iptr_row++)>=0)
	{
	  while ((*iptr_list>=0)&&(elm>*iptr_list))
	    {j++; iptr_list++;}
	  if (*iptr_list<0) break;
	  if (elm == *iptr_list)
	    {
	      *iptr_ja++=j++;
	      if (elm==k)
		{v_hold = dptr;}
	      else
		{
		  *dptr = *dptr_vals;
		  row_sum += *dptr;
		}

	      /* old */
	      /* = (double) (*(iptr_ia+1) - *(iptr_ia) - 1); iptr_ia++;}*/
	      /* for spd */
	      /*
	      {*dptr = (double) (*(iptr_ia+1) - *(iptr_ia)); iptr_ia++;}
	      */

	      dptr++;
	      iptr_list++;
	    }
	  dptr_vals++;
	}
      if (row_sum==0.0)
	{
	  for (i=0; i<n; i++)
	    {pre_nek_color[list[i]] = 15;}
	  /* drmesh_color_(); */
	  pause();
	  error_msg_warning("lap_sub_matrix() :: definitely disconnected!\n");
	}

      *v_hold = -row_sum;
      wt+= -row_sum;
    }

  check_row_sum(n,*ia,*ja,*vals);

#ifdef DEBUG
  error_msg_warning("lap_sub_matrix() :: end\n");
#endif

  return(wt);

  /* later :: chk to make sure graph is connected */

#ifdef DEBUG_0
  iptr_ia = ia+1;
  iptr_ja = ja;
  dptr   = vals;
  for (i=0; i<n; iptr_ia++, i++)
    {
      k = *(iptr_ia) - *(iptr_ia-1);
      for (j=0; j<k; j++)
	{printf("%3d ",*iptr_ja++);}
      printf("\n");

      for (j=0; j<k; j++)
	{printf("%3.0f ",*dptr++);}
      printf("\n");
    }
#endif
}



/********************************sparse_matrix.c*******************************
Function: 

Input : 
Output: 
Return: 
Description:  

Make sure that row sum of laplacian matrix is in fact zero as it should be.

*********************************sparse_matrix.c******************************/
void
check_row_sum(int n, int *ia, int *ja, double *vals)
{
  int i;
  double *b, *x;


#ifdef DEBUG
  error_msg_warning("check_row_sum() :: begin\n");
#endif

  b = (double *) bss_malloc(n*sizeof(double));
  x = (double *) bss_malloc(n*sizeof(double));

  dvec_set(x,1.0,n);
  matrix_mult(b,x,n,ia,ja,vals);

  for (i=0; i<n; i++)
    {
      if (fabs(b[i])>1.0e-14) 
	{error_msg_fatal("check_row_sum() :: row sum %d = %g\n",i);}
    }


  bss_free(b);
  bss_free(x);

#ifdef DEBUG
  error_msg_warning("check_row_sum() :: end\n");
#endif
}



/********************************sparse_matrix.c*******************************
Function: restricted_kl()

Input : 
Output: 
Return: 
Description:  

   well you know ... I'll get to it!
*********************************sparse_matrix.c******************************/
void
restricted_kl(int *list, int n, int *ia, int *ja, double *vals, int *color)
{
  ;
}



/********************************sparse_matrix.c*******************************
Function: brute_separate_graph()

Input : 
Output: 
Return: 
Description:  

bisect graph and highlight those which have elements with cut edges
   o <0 lhs
   o -1 lhs element has edge cut
   o >0 rhs
   o +1 lhs element has edge cut

   ***note we're assuming that the graphs are connected***

*********************************sparse_matrix.c******************************/
int *
brute_sep_graph(int *list, int n, int *ia, int *ja, double *vals, double *rhs, 
		int *ns)
{
  int i,j,k;
  int *color;
  double *x, *b, rs, nc=0.0;
  double *fiedler;
  int *map;
  int nnz, uh_nnz, bw;
  int *iptr, sm_in[32], sm_out[32];
  int v1,v2,v3,v4;
  double d12,d13,d14,d23,d24,d34;
  int who;
  int sv[6];
  double dist[6];
  int passed = TRUE;


#ifdef DEBUG
  error_msg_warning("brute_sep_graph() :: begin\n");
#endif

  color = (int *) bss_malloc(n*INT_LEN);

  switch (n) {
  case 1:
    color[0] = -1; return(color);
    break;
  case 2:
    color[0] = -1; color[1] =  1; return(color);
    break;
  case 3:
    /* fully connected so choose closest pair */
    /* if (ia[n]==(n*n)) */
    if (0)
      {
	v1 = list[0];
	v2 = list[1];
	v3 = list[2];

	ivec_c_index(sv,3);
	if (rsb_dim==_2D)
	  {
	    d12  = pow(xc[v1]-xc[v2],2.0);
	    d12 += pow(yc[v1]-yc[v2],2.0);
	    dist[0] = d12 = sqrt(d12);
	    d13  = pow(xc[v1]-xc[v3],2.0);
	    d13 += pow(yc[v1]-yc[v3],2.0);
	    dist[1] = d13 = sqrt(d13);
	    d23  = pow(xc[v2]-xc[v3],2.0);
	    d23 += pow(yc[v2]-yc[v3],2.0);
	    dist[2] = d23 = sqrt(d23);
	  }
	else
	  {
	    d12  = pow(xc[v1]-xc[v2],2.0);
	    d12 += pow(yc[v1]-yc[v2],2.0);
	    d12 += pow(zc[v1]-zc[v2],2.0);
	    dist[0] = d12 = sqrt(d12);
	    d13  = pow(xc[v1]-xc[v3],2.0);
	    d13 += pow(yc[v1]-yc[v3],2.0);
	    d13 += pow(zc[v1]-zc[v3],2.0);
	    dist[1] = d13 = sqrt(d13);
	    d23  = pow(xc[v2]-xc[v3],2.0);
	    d23 += pow(yc[v2]-yc[v3],2.0);
	    d23 += pow(zc[v2]-zc[v3],2.0);
	    dist[2] = d23 = sqrt(d23);
	  }

	dvec_sort_companion(dist,sv,3);
	switch (sv[0]) {
	case 0:
	  {
	    color[0] = -1;
	    color[1] = -1;
	    color[2] =  1;
	    break;
	  }
	case 1:
	  {
	    color[0] = -1;
	    color[1] =  1;
	    color[2] = -1;

	    break;
	  }
	case 2:
	  {
	    color[0] =  1;
	    color[1] = -1;
	    color[2] = -1;
	    break;
	  }
	default:
	  {
	    error_msg_warning("brute_sep_graph() :: bad order n=3\n");
	    passed=FALSE;
	    break;
	  }
	}
      }
    /* not fully connected so use mindegree ordering to split */
    else
      {
	nnz = ia[n];
	uh_nnz = nnz - n;
	if (uh_nnz&1)
	  {
	    error_msg_warning("brute_sep_graph_() :: matrix not symmetric!!\n");
	    passed=FALSE;
	    break;
	  }

	if (uh_nnz>0)
	  {uh_nnz>>=1;}
	iptr = sm_in;
	*iptr++ = n;
	*iptr++ = k = uh_nnz;
	for (i=j=0; i<n; i++)
	  {
	    while (j<ia[i+1])
	      {
		if (ja[j]>i)
		  {*iptr++ = i+1; *iptr++ = ja[j]+1; k--;}
		j++;
	      }
	  }
	if (k)
	  {
	    error_msg_fatal("brute_sep_graph () :: pair list failure!\n");
	    passed=FALSE;
	    break;
	  }

#ifdef DEBUG	  
	error_msg_warning("A LIST :: ");
	for (i=0; i<n; i++)
	  {printf("%d ",list[i]);}
	printf("\n");
	fflush(stdout);
	dump_laplacian(n,ia,ja,vals);
	error_msg_warning("done dump_laplacian\n");
#endif

	/* use min degree ordering */
	bw=sm_bandwidth_reduction(sm_in,sm_out,0);
	error_msg_warning("A LIST bw = %d\n",bw);

	/* oops */
	if (bw<=0)
	  {
	    error_msg_warning("B LIST :: ");
	    for (i=0; i<n; i++)
	      {printf("%d ",list[i]);}
	    printf("\n");
	    fflush(stdout);
	    dump_laplacian(n,ia,ja,vals);
	    pre_nek_color[0] = 15;
	    pre_nek_color[1] = 15;
	    pre_nek_color[2] = 15;
	    /* drmesh_color_(); */
	    error_msg_warning("brute_sep_graph_ () :: sm failured for n=%d!\n",n);
	    passed=FALSE;
	    break;
	  }

	color[sm_out[0]-1] = -1;
	color[sm_out[1]-1] =  1;
	color[sm_out[2]-1] =  1;
      }
    break;
  case 4:
    /* fully connected so choose closest pair */
    /* if (ia[n]==(n*n)) */
    if (0)
      {
#ifdef DEBUG
	error_msg_warning("fully connected!\n");
#endif
	v1 = list[0];
	v2 = list[1];
	v3 = list[2];
	v4 = list[3];

	ivec_c_index(sv,6);

	if (rsb_dim==_2D)
	  {
	    d12  = pow(xc[v1]-xc[v2],2.0);
	    d12 += pow(yc[v1]-yc[v2],2.0);
	    dist[0] = d12 = sqrt(d12);
	    d13  = pow(xc[v1]-xc[v3],2.0);
	    d13 += pow(yc[v1]-yc[v3],2.0);
	    dist[1] = d13 = sqrt(d13);
	    d14  = pow(xc[v1]-xc[v4],2.0);
	    d14 += pow(yc[v1]-yc[v4],2.0);
	    dist[2] = d14 = sqrt(d14);
	    d23  = pow(xc[v2]-xc[v3],2.0);
	    d23 += pow(yc[v2]-yc[v3],2.0);
	    dist[3] = d23 = sqrt(d23);
	    d24  = pow(xc[v2]-xc[v4],2.0);
	    d24 += pow(yc[v2]-yc[v4],2.0);
	    dist[4] = d24 = sqrt(d24);
	    d34  = pow(xc[v3]-xc[v4],2.0);
	    d34 += pow(yc[v3]-yc[v4],2.0);
	    dist[5] = d34 = sqrt(d34);
	  }
	else
	  {
	    d12  = pow(xc[v1]-xc[v2],2.0);
	    d12 += pow(yc[v1]-yc[v2],2.0);
	    d12 += pow(zc[v1]-zc[v2],2.0);
	    dist[0] = d12 = sqrt(d12);
	    d13  = pow(xc[v1]-xc[v3],2.0);
	    d13 += pow(yc[v1]-yc[v3],2.0);
	    d13 += pow(zc[v1]-zc[v3],2.0);
	    dist[1] = d13 = sqrt(d13);
	    d14  = pow(xc[v1]-xc[v4],2.0);
	    d14 += pow(yc[v1]-yc[v4],2.0);
	    d14 += pow(zc[v1]-zc[v4],2.0);
	    dist[2] = d14 = sqrt(d14);
	    d23  = pow(xc[v2]-xc[v3],2.0);
	    d23 += pow(yc[v2]-yc[v3],2.0);
	    d23 += pow(zc[v2]-zc[v3],2.0);
	    dist[3] = d23 = sqrt(d23);
	    d24  = pow(xc[v2]-xc[v4],2.0);
	    d24 += pow(yc[v2]-yc[v4],2.0);
	    d24 += pow(zc[v2]-zc[v4],2.0);
	    dist[4] = d24 = sqrt(d24);
	    d34  = pow(xc[v3]-xc[v4],2.0);
	    d34 += pow(yc[v3]-yc[v4],2.0);
	    d34 += pow(zc[v3]-zc[v4],2.0);
	    dist[5] = d34 = sqrt(d34);
	  }
	  
	dvec_sort_companion(dist,sv,6);
	switch (sv[0]) {
	case 0:
	  {
	    color[0] = -1;
	    color[1] = -1;
	    color[2] =  1;
	    color[3] =  1;
	    break;
	  }
	case 1:
	  {
	    color[0] = -1;
	    color[1] =  1;
	    color[2] = -1;
	    color[3] =  1;
	    break;
	  }
	case 2:
	  {
	    color[0] = -1;
	    color[1] =  1;
	    color[2] =  1;
	    color[3] = -1;
	    break;
	  }
	case 3:
	  {
	    color[0] =  1;
	    color[1] = -1;
	    color[2] = -1;
	    color[3] =  1;
	    break;
	  }
	case 4:
	  {
	    color[0] =  1;
	    color[1] = -1;
	    color[2] =  1;
	    color[3] = -1;
	    break;
	  }
	case 5:
	  {
	    color[0] =  1;
	    color[1] =  1;
	    color[2] = -1;
	    color[3] = -1;
	    break;
	  }
	default:
	  {
	    error_msg_warning("brute_sep_graph() :: bad order n=4\n");
	    passed=FALSE;
	    break;
	  }
	}
      }
    else
      {
#ifdef DEBUG
	error_msg_warning("not fully connected!\n");
#endif
	nnz = ia[n];
	uh_nnz = nnz - n;
	if (uh_nnz&1)
	  {
	    error_msg_fatal("brute_sep_graph_() :: matrix not symmetric!!\n");
	    passed=FALSE;
	    break;
	  }
	if (uh_nnz>0)
	  {uh_nnz>>=1;}
	ivec_zero(sm_in,32);
	ivec_zero(sm_out,32);
	iptr = sm_in;
	*iptr++ = n;
	*iptr++ = k = uh_nnz;
#ifdef DEBUG
	error_msg_warning("%d, %d :: ",n,k);
#endif
	for (i=j=0; i<n; i++)
	  {
	    while (j<ia[i+1])
	      {
		if (ja[j]>i)
		  {
		    *iptr++ = i+1; *iptr++ = ja[j]+1; k--;
#ifdef DEBUG
		    error_msg_warning("(%d,%d) ",i+1,ja[j]+1);
#endif
		  }
		j++;
	      }
	  }
#ifdef DEBUG
	printf("\n");
#endif

	if (k)
	  {
	    error_msg_warning("brute_sep_graph () :: pair list failure!\n");
	    printf("k=%d\n",k);
	    printf("C LIST :: ");
	    for (i=0; i<n; i++)
	      {printf("%d ",list[i]);}
	    printf("\n");
	    fflush(stdout);
	    dump_laplacian(n,ia,ja,vals);
	    passed=FALSE;
	    break;
	  }

#ifdef DEBUG
	error_msg_warning("D LIST :: ");
	for (i=0; i<n; i++)
	  {printf("%d ",list[i]);}
	printf("\n");
	fflush(stdout);
	dump_laplacian(n,ia,ja,vals);
#endif

#ifdef DEBUG
	error_msg_warning("calling sm_...() ... ");
#endif
	/* use min degree ordering */
	bw=sm_bandwidth_reduction(sm_in,sm_out,0);

#ifdef DEBUG
	error_msg_warning("... returning from sm_...()!\n");
#endif
	  
	if (bw<=0)
	  {
	    error_msg_warning("E LIST :: ");
	    for (i=0; i<n; i++)
	      {printf("%d ",list[i]);}
	    printf("\n");
	    fflush(stdout);
	    dump_laplacian(n,ia,ja,vals);
	    pre_nek_color[list[0]] = 15;
	    pre_nek_color[list[1]] = 15;
	    pre_nek_color[list[2]] = 15;
	    pre_nek_color[list[3]] = 15;
	    /* drmesh_color_(); */
	    error_msg_warning("brute_sep_graph_ () :: sm_...() failure!\n");
	    pause();
	    passed=FALSE;
	    break;
	  }

#ifdef DEBUG
	error_msg_warning("OUT :: ");
	for (i=0; i<n; i++)
	  {printf("%d ",sm_out[i]-1);}
	printf("\n");
#endif
	v1 = sm_out[0]-1;
	v2 = sm_out[1]-1;
	v3 = sm_out[2]-1;
	v4 = sm_out[3]-1;


	/* make sure v1 and v2 are connected */
	k=FALSE;
	for (i=ia[v1]; i<ia[v1+1]; i++)
	  { 
	    if (ja[i]==v2)
	      {k=TRUE; break;}
	  }

#ifdef DEBUG
	if (!k)
	  {
	    error_msg_warning("F LIST :: ");
	    for (i=0; i<n; i++)
	      {printf("%d ",list[i]);}
	    printf("\n");
	    fflush(stdout);
	    dump_laplacian(n,ia,ja,vals);
	    pre_nek_color[list[v1]] = 2;
	    pre_nek_color[list[v2]] = 2;
	    /* drmesh_color_(); */
	    pause();
	  }

	/* something is not right! */
	color[v1] = -1;
	color[v2] = -1;
	color[v3] =  1;
	color[v4] =  1;
	break;
#endif

	/* want v1,v2 unless v1 connected to v3 and d(v1,v3) < d(v1,v2) */
	k=FALSE;
	for (i=ia[v1]; i<ia[v1+1]; i++)
	  { 
	    if (ja[i]==v3)
	      {k=TRUE; break;}
	  }

	/* v1 only connected to v2 */
	/* if (!k)*/
	if (1)
	  {
	    color[v1] = -1;
	    color[v2] = -1;
	    color[v3] =  1;
	    color[v4] =  1;
	  }
	/* v1 connected to v3 - note that it can't be connected to v4 */
	else
	  {
	    if (rsb_dim==_2D)
	      {
		d12  = pow(xc[list[v1]]-xc[list[v2]],2.0);
		d12 += pow(yc[list[v1]]-yc[list[v2]],2.0);
		d12 = sqrt(d12);
		d13  = pow(xc[list[v1]]-xc[list[v3]],2.0);
		d13 += pow(yc[list[v1]]-yc[list[v3]],2.0);
		d13 = sqrt(d13);
	      }
	    else
	      {
		d12  = pow(xc[list[v1]]-xc[list[v2]],2.0);
		d12 += pow(yc[list[v1]]-yc[list[v2]],2.0);
		d12 += pow(zc[list[v1]]-zc[list[v2]],2.0);
		d12 = sqrt(d12);
		d13  = pow(xc[list[v1]]-xc[list[v3]],2.0);
		d13 += pow(yc[list[v1]]-yc[list[v3]],2.0);
		d13 += pow(zc[list[v1]]-zc[list[v3]],2.0);
		d13 = sqrt(d13);
	      }

	    if (d13<d12)
	      {
		color[v1] = -1;
		color[v2] =  1;
		color[v3] = -1;
		color[v4] =  1;
	      }
	    else
	      { 
 		k=FALSE;
		color[v1] = -1;
		color[v2] = -1;
		color[v3] =  1;
		color[v4] =  1;
	      }
	  }


	/* k==true  ==> (v1,v3) so check (v2,v4) */
	if (0)
	  {
	    passed=FALSE;
	    for (i=ia[v2]; i<ia[v2+1]; i++)
	      { 
		if (ja[i]==v4)
		  {passed=TRUE; break;}
	      }
	  }
	/* k==false ==> (v1,v2) so check (v3,v4) */
	else
	  {
	    passed=FALSE;
	    for (i=ia[v3]; i<ia[v3+1]; i++)
	      { 
		if (ja[i]==v4)
		  {passed=TRUE; break;}
	      }
	  }
      }
    break;
  default:
    error_msg_warning("brute_sep_graph() :: can't handle %d elements!\n",n);
    passed=FALSE;
    break;
  }

  if (passed==FALSE)
    {bss_free(color); color=NULL;}

#ifdef DEBUG
  error_msg_warning("brute_sep_graph() :: end\n");
  error_msg_warning("passed=%d\n",passed);
#endif

  return(color);
}


/********************************sparse_matrix.c*******************************
Function: separate_graph()

Input : 
Output: 
Return: 
Description:  

bisect graph and highlight those which have elements with cut edges
   o <0 lhs
   o -2 lhs element has edge cut
   o >0 rhs
   o +2 lhs element has edge cut


*********************************sparse_matrix.c******************************/
int *
sep_graph(int *list, int n, int *ia, int *ja, double *vals, double *rhs, int *ns)
{
  int i,j,k;
  int *color;
  double *x, *b, rs, nc=0.0;
  double *fiedler;
  int *map, *sm_in, *sm_out, *iptr;
  int nprs, nnz, nl, nr, rhs_bw, lhs_bw;
  int row, col;
  int connected, passes;


#ifdef DEBUG
  error_msg_warning("sep_graph() :: begin\n");
#endif

  /* brute separate */
  if (n<=4)
    {
      if (color=brute_sep_graph(list,n,ia,ja,vals,rhs,ns))
	{
	  x = (double *) bss_malloc(n*sizeof(double));
	  b = (double *) bss_malloc(n*sizeof(double));
	  dvec_zero(x,n);

	  for (i=j=0; i<n; i++)
	    {
	      if (color[i]>0)
		{x[i]=1.0; j++;}
	    }

	  if (abs(n/2-j)>1)
	    {
	      printf("COLOR :: ");
	      for (i=0; i<n; i++)
		{printf("%d ",color[i]);}
	      printf("\nG LIST :: ");
	      for (i=0; i<n; i++)
		{printf("%d ",list[i]);}
	      printf("\n");
	      fflush(stdout);
	      dump_laplacian(n,ia,ja,vals);
	      fflush(stdout);
	      error_msg_fatal("sep_graph() :: brute failed!\n");
	    }
      
	  matrix_mult(b,x,n,ia,ja,vals);
	  for (i=0; i<n; i++)
	    {
	      if ((rs=fabs(b[i]))>0.0)
		{color[i]*=2; nc += rs;}
	    }
	  nc+=0.5;
	  *ns = (int) nc;

	  bss_free(x);
	  bss_free(b);
	}
      else
	{*ns=INT_MAX;}
      return(color);
    }
    
  /* general case */
  i = ia[n]-n;
  if (i&1) 
    {
      printf("\nH LIST :: ");
      for (i=0; i<n; i++)
	{printf("%d ",list[i]);}
      printf("\n");
      fflush(stdout);
      dump_laplacian(n,ia,ja,vals);
      fflush(stdout);
      error_msg_fatal("sep_graph_() :: matrix not symmetric!!\n");
    }

  sm_in  = iptr = (int *) bss_malloc((i+2)*INT_LEN);
  sm_out = (int *) bss_malloc(n*INT_LEN);
  color = (int *) bss_malloc(n*INT_LEN);
  map   = (int *) bss_malloc(n*INT_LEN);
  x = (double *) bss_malloc(n*sizeof(double));
  b = (double *) bss_malloc(n*sizeof(double));
  connected = TRUE;

  nl=n/2;
  nr=n-nl;
  lhs_bw=rhs_bw=0;
  passes=1;
  for(;;)
    {
      /* get fiedler vector associated w/subgraph */
      dvec_dump(rhs,n,passes,nl,"rhs");
      if (fiedler=lanczos(rhs, n, ia, ja, vals))
	{
	  dvec_dump(fiedler,n,passes,nl,"fiedler");
	  ivec_c_index(map,n);
	  ivec_zero(color,n);
	  dvec_zero(x,n);

	  /* determine lhs,rhs */
	  dvec_sort_companion(fiedler,map,n);
	  bss_free(fiedler);
	  for (i=0; i<n/2; i++)
	    {color[map[i]] = -1; }

	  /* process lhs */
	  *iptr++ = nl;
	  *iptr++ = k = 0;
	  
	  /* I should build a new map and indirect but this will have to wait */
	  for (i=j=0; i<n; i++)
	    {
	      row = ivec_linear_search(i,map,nl);
	      if (row!=-1)
		{
		  row++;
		  while (j<ia[i+1])
		    {
		      col = ivec_linear_search(ja[j],map,nl);
		      if (col!=-1)
			{
			  col++;
			  printf("(%d,%d) ",row,col);
			  if (row<col)
			    {
			      *iptr++ = col; 
			      *iptr++ = row;
			      k++;
#ifdef DEBUG			      
			      if (n<=32)
				{printf("(%d,%d) ",row,col);}
#endif
			    }
			}
		      j++;
		    }
#ifdef DEBUG			      
		  if (n<=32)
		    {printf("\n");}
#endif
		}
	      else
		{j=ia[i+1];}
	    }
	  sm_in[1] = k;
	  
#ifdef DEBUG			      
	  if (!k)
	    {error_msg_fatal("sep_graph () :: no pairs to lhs?\n");}
#endif

	  /* check to make sure that potential lhs sub-graph is connected */
	  /* use elm 1 to orderng ... we're throwing it away anyway */
	  lhs_bw=sm_bandwidth_reduction(sm_in,sm_out,1);

#ifdef DEBUG			      
	  if (n<=32)
	    {printf("sep_graph () :: nl=%d, nprs=%d, bw=%d\n",nl,k,lhs_bw);}
#endif

	  /* process rhs */
	  for (i=n/2; i<n; i++)
	    {color[map[i]] =  1; x[map[i]]=1.0;}

	  iptr = sm_in;	  
	  *iptr++ = nr;
	  *iptr++ = k = 0;	      
	  /* I should build a new map and indirect but this will have to wait */
	  for (i=j=0; i<n; i++)
	    {
	      row = ivec_linear_search(i,map+nl,nr);
	      if (row!=-1)
		{
		  row++;
		  while (j<ia[i+1])
		    {
		      col = ivec_linear_search(ja[j],map+nl,nr);
		      if (col!=-1)
			{
			  col++;
			  if (row<col)
			    {
			      *iptr++ = col; 
			      *iptr++ = row;
			      k++;
#ifdef DEBUG			      
			      if (n<=32)
				{printf("(%d,%d) ",row,col);}
#endif
			    }
			}
		      j++;
		    }
#ifdef DEBUG			      
		  if (n<=32)
		    {printf("\n");}
#endif
		}
	      else
		{j=ia[i+1];}
	    }
	  sm_in[1] = k;	      

#ifdef DEBUG			      
	  if (!k)
	    {error_msg_fatal("sep_graph () :: no pairs to rhs?\n");}
#endif

	  /* check to make sure potnetial rhs sub-graph is connected */
	  /* use elm 1 to orderng ... we're throwing it away anyway */
	  rhs_bw=sm_bandwidth_reduction(sm_in,sm_out,1);
	
#ifdef DEBUG			      
	  if (n<=32)
	    {printf("sep_graph () :: nr=%d, nprs=%d, bw=%d\n",nr,k,lhs_bw);}
#endif
	}
      else
	{
#ifdef DEBUG
	  error_msg_warning("sep_graph() :: lanczos failure!\n");
#endif
	}

      /* connected ==> done! */
      if ((rhs_bw>0)&&(lhs_bw>0))
	{
#ifdef DEBUG
	  error_msg_warning("lhs_bw=%d, rhs_bw=%d\n",lhs_bw,rhs_bw); 
#endif
	  break;
	}
	  
      /* I'm not willing to wait forever!                  */
      /* actually ... this isn't working properly! 6.19.99 */
      if (passes>0)
	{

#ifdef DEBUG
	  error_msg_warning("sep_graph() :: too many passes! (%d,%d)\n",lhs_bw,rhs_bw); 
	  printf("I LIST :: ");
	  for (i=0; i<n; i++)
	    {printf("%d ",list[i]);}
	  printf("\n");
	  fflush(stdout);
	  dump_laplacian(n,ia,ja,vals);
	  fflush(stdout);

	  /*
	  refresh_();
	  drgrid_();
	  */

	  k = 2;
	  for (i=0; i<n; i++)
	    {
	      j = list[i]+1;
	      /* hmt_drawel_color_(&j,&k); */
	    }
	  fflush(stdout);
#endif
	  connected=FALSE;
	  break; 
	}

      /* not connected but we'll try again w/diff rhs */
      if (passes%5)
	{
#ifdef DEBUG
	  error_msg_warning("sep_graph() :: perturbing rhs!\n"); 
#endif
	  hmt_perturb(rhs,n);
	}
      else
	{
#ifdef DEBUG
	  error_msg_warning("sep_graph() :: randomizing rhs!\n"); 
#endif
	  Init_rhs_random(rhs,n);
	}

      passes++;
    }
#ifdef DEBUG      
  error_msg_warning("sep_graph () :: num passes=%d\n",passes);
#endif
  if (connected)
    {
      error_msg_warning("sep_graph () :: connected\n");
      /* local clean up to minimize number of cut edges */
      restricted_kl(list,n,ia,ja,vals,color);

      /* mark elms w/edges cut */
      /* need to fix x if restricted_kl does something!!! */
      matrix_mult(b,x,n,ia,ja,vals);
      nc = 0.0;
      for (i=0; i<n; i++)
	{
	  if ((rs=fabs(b[i]))>0.0)
	    {color[i]*=2; nc += rs;}
	}
	  
      nc+=0.5;
      *ns = (int) nc;	

      ivec_dump(color,n,*ns,passes,"color in sep_graph");
    }
  else
    {
      error_msg_warning("sep_graph () :: not connected\n");
      *ns = INT_MAX;
      bss_free(color);
      color=NULL;
    }

  /* return scratch space */
  bss_free(x);
  bss_free(b);
  bss_free(map);
  bss_free(sm_in);
  bss_free(sm_out);

#ifdef DEBUG
  error_msg_warning("sep_graph() :: end\n");
#endif

  return(color);
}



/********************************sparse_matrix.c*******************************
Function: label_outflow()

Input : 
Output: 
Return: 
Description:  
*********************************sparse_matrix.c******************************/
int *
label_outflow(int *vertex, int nv, int *num_unique, int *num_outflow)
{
  int i, j, off, nvo, nvu;
  int *iptr;


#ifdef DEBUG
  error_msg_warning("label_outflow() :: begin\n");
#endif

  *num_unique  = nvu = vertex[nv];
  *num_outflow = nvo = vertex[nv+1];

#ifdef DEBUG
  printf("nv=%d, nu=%d, no=%d\n",nv,nvu,nvo);
#endif

  iptr = (int *)bss_malloc(INT_LEN*(nvu+1));
  ivec_zero(iptr,nvu+1);

  for (i=j=0; i<nv; i++)
    {
      if (vertex[i] < 0)
	{
	  off = vertex[i] = -vertex[i];
	  if (!iptr[off])
	    {
	      iptr[off] = nvu-j;
	      j++;
	    }
	}
    }

  if (j!=nvo)
    {error_msg_fatal("label_outflow() :: found=%d, actual=%d\n",j,nvo);}

#ifdef DEBUG
  error_msg_warning("label_outflow() :: end\n");
#endif

  return(iptr);
}



/********************************sparse_matrix.c*******************************
Function: det_sep()

Input : 
Output: 
Return: 
Description:  
*********************************sparse_matrix.c******************************/
int *
det_sep(int *list, int *color, int n, int *vertex, int nv, int *vertex_order, 
	int nvu, int ub, int *num_sep)
{
  int i, j, nc, off, ns=0, nsm, v;
  int *iptr=NULL;

#ifdef DEBUG
  error_msg_warning("det_sep() :: begin\n");
#endif
  
  /* nsm should be an upper bound on the number of vertices to be labelled */
  nsm = n;
  nc = (rsb_dim==_2D) ? 4 : 8;
  nsm *= nc;
  nsm = MAX(nsm,nc);

  iptr = (int *) bss_malloc(INT_LEN*nsm);

  if (color==NULL)
    {
      for (ns=i=0; i<n; i++)
	{
	  off = nc*list[i];
	  for (j=0; j<nc; j++)
	    {
	      v = vertex[off+j];
	      if (vertex_order[v] == 0)
		{
#ifdef DEBUG		  		  
		  printf("elm=%d,v=%d\n",off,v);
#endif
		  if (ns>=nsm)
		    {error_msg_fatal("oops ns about to exceed nsm=%d!\n",nsm);}
		  iptr[ns] = vertex_order[v] = ub-ns;
		  ns++;
		}
	    }
	}
    }
  else
    {
      for (i=0; i<n; i++)
	{
	  if (color[i] == -2)
	    {
	      off = nc*list[i];
	      for (j=0; j<nc; j++)
		{
		  v = vertex[off+j];
		  if (vertex_order[v] <= 0)
		    {
		      /* printf("elm=%d,v=%d\n",off,v); */
		      vertex_order[v]--;
		    }
		}
	    }
	}


      for (ns=i=0; i<n; i++)
	{
	  if (color[i] == 2)
	    {
	      off = nc*list[i];
	      for (j=0; j<nc; j++)
		{
		  v = vertex[off+j];
		  if (vertex_order[v] < 0)
		    {
		      /* printf("elm=%d,v=%d\n",off,v); */
		      if (ns>=nsm)
			{error_msg_fatal("oops ns about to exceed nsm=%d!\n",nsm);}

		      iptr[ns] = vertex_order[v] = ub-ns;
		      ns++;
		    }
		}
	    }
	}
    }

  /* reset the unlabelled vertices to 0 state */
  for (i=1; i<=nvu; i++)
    {if (vertex_order[i]<0) {vertex_order[i]=0;}}

  *num_sep = ns;
  ivec_sort(iptr,ns);

#ifdef DEBUG
  error_msg_warning("det_sep() :: end\n");
#endif

  return(iptr);
}



/********************************sparse_matrix.c*******************************
Function: det_lc()

Input : 
Output: 
Return: 
Description:  
*********************************sparse_matrix.c******************************/
struct tree_node *
det_lc_vertex(int *list, int *color, double *rhs, int n, int *vertex, int nv, 
       int *vertex_order, int nvu, int lb, int *num_verts)
{
  int i, j, k;
  int nl,nc,v,off;
  int *iptr;
  double *dptr;
  struct tree_node *t_ptr;
  
#ifdef DEBUG
  error_msg_warning("det_lc() :: begin\n");
#endif

  t_ptr = new_tree_node();
  t_ptr->active = NULL;
  t_ptr->rhs    = NULL;
  nc = (rsb_dim==_2D) ? 4 : 8;
  for (nl=i=0; i<n; i++)
    {
      if (color[i] < 0)
	{ 
	  nl++;
	  off = nc*list[i];
	  for (j=0; j<nc; j++)
	    {
	      v = vertex[off+j];
	      if (vertex_order[v] <= 0)
		{vertex_order[v]--;}
	    }
	}
    }

  if (!nl)
    {printf("det_lc() :: nl=%d\n",nl);}
  t_ptr->num_act = nl;
  t_ptr->active = iptr = (int *) bss_malloc((nl+1)*INT_LEN);
  t_ptr->rhs    = dptr = (double *) bss_malloc((nl+1)*sizeof(double));
  for (nl=i=0; i<n; i++)
    {
      if (color[i] < 0)
	{
	  nl++; 
	  *iptr++ = list[i];
	  *dptr++ = rhs[i];
	}
    }
  *iptr = -1;

  if (nl!=t_ptr->num_act)
    {error_msg_fatal("lc ct off!");}

  /* reset the unlabelled vertices to 0 state */
  for (nc=0,i=1; i<=nvu; i++)
    {if (vertex_order[i]<0) {nc++; vertex_order[i]=0;}}


  t_ptr->lb = lb;
  *num_verts = nc;
  t_ptr->ub  = nc+lb-1;

#ifdef DEBUG
  error_msg_warning("det_lc() :: end\n");
#endif

  return(t_ptr);
}


/********************************sparse_matrix.c*******************************
Function: det_lc()

Input : 
Output: 
Return: 
Description:  
*********************************sparse_matrix.c******************************/
struct tree_node *
det_lc(int *list, int *color, double *rhs, int n)
{
  int i, j, k;
  int nl,nc,v,off;
  int *iptr;
  double *dptr;
  struct tree_node *t_ptr;
  
#ifdef DEBUG
  error_msg_warning("det_lc() :: begin\n");
#endif

  t_ptr = new_tree_node();
  t_ptr->active = NULL;
  t_ptr->rhs    = NULL;

  for (nl=i=0; i<n; i++)
    {
      if (color[i] < 0)
	{nl++;}
    }

  if (!nl)
    {printf("det_lc() :: nl=%d\n",nl);}

  t_ptr->num_act = nl;
  t_ptr->active = iptr = (int *) bss_malloc((nl+1)*INT_LEN);
  t_ptr->rhs    = dptr = (double *) bss_malloc((nl+1)*sizeof(double));
  for (nl=i=0; i<n; i++)
    {
      if (color[i] < 0)
	{
	  nl++; 
	  *iptr++ = list[i];
	  *dptr++ = rhs[i];
	}
    }
  *iptr = -1;

  if (nl!=t_ptr->num_act)
    {error_msg_fatal("lc ct off!");}

#ifdef DEBUG
  error_msg_warning("det_lc() :: end w/%d of %d going left\n",nl,n);
#endif

  return(t_ptr);
}


/********************************sparse_matrix.c*******************************
Function: det_rc()

Input : 
Output: 
Return: 
Description:  
*********************************sparse_matrix.c******************************/
struct tree_node *
det_rc_vertex(int *list, int *color, double *rhs, int n, int *vertex, int nv, 
       int *vertex_order, int nvu, int lb, int *num_verts)
{
  int i, j, k;
  int nl,nc,v,off;
  int *iptr;
  double *dptr;
  struct tree_node *t_ptr;
  
#ifdef DEBUG
  error_msg_warning("det_rc() :: begin\n");
#endif

  t_ptr = new_tree_node();
  t_ptr->active = NULL;
  t_ptr->rhs    = NULL;
  nc = (rsb_dim==_2D) ? 4 : 8;
  for (nl=i=0; i<n; i++)
    {
      if (color[i] > 0)
	{
	  nl++;
	  off = nc*list[i];
	  for (j=0; j<nc; j++)
	    {
	      v = vertex[off+j];
	      if (vertex_order[v] <= 0)
		{vertex_order[v]--;}
	    }
	}
    }

  if (!nl)
    {printf("det_rc() :: nl=%d\n",nl);}
  t_ptr->num_act = nl;
  t_ptr->active = iptr = (int *) bss_malloc((nl+1)*INT_LEN);
  t_ptr->rhs    = dptr = (double *) bss_malloc((nl+1)*sizeof(double));
  for (nl=i=0; i<n; i++)
    {
      if (color[i] > 0)
	{
	  nl++; 
	  *iptr++ = list[i];
	  *dptr++ = rhs[i];
	}
    }
  *iptr = -1;

  if (nl!=t_ptr->num_act)
    {error_msg_fatal("lc ct off!");}

  /* reset the unlabelled vertices to 0 state */
  for (nc=0,i=1; i<=nvu; i++)
    {if (vertex_order[i]<0) {nc++; vertex_order[i]=0;}}


  t_ptr->lb = lb;
  *num_verts = nc;
  t_ptr->ub  = nc+lb-1;

#ifdef DEBUG
  error_msg_warning("det_rc() :: end\n");
#endif

  return(t_ptr);
}


/********************************sparse_matrix.c*******************************
Function: det_rc()

Input : 
Output: 
Return: 
Description:  
*********************************sparse_matrix.c******************************/
struct tree_node *
det_rc(int *list, int *color, double *rhs, int n)
{
  int i, j, k;
  int nl,nc,v,off;
  int *iptr;
  double *dptr;
  struct tree_node *t_ptr;
  
#ifdef DEBUG
  error_msg_warning("det_rc() :: begin\n");
#endif

  t_ptr = new_tree_node();
  t_ptr->active = NULL;
  t_ptr->rhs    = NULL;

  for (nl=i=0; i<n; i++)
    {
      if (color[i] > 0)
	{nl++;}
    }

  if (!nl)
    {printf("det_rc() :: nl=%d\n",nl);}
  t_ptr->num_act = nl;
  t_ptr->active = iptr = (int *) bss_malloc((nl+1)*INT_LEN);
  t_ptr->rhs    = dptr = (double *) bss_malloc((nl+1)*sizeof(double));
  for (nl=i=0; i<n; i++)
    {
      if (color[i] > 0)
	{
	  nl++; 
	  *iptr++ = list[i];
	  *dptr++ = rhs[i];
	}
    }
  *iptr = -1;

  if (nl!=t_ptr->num_act)
    {error_msg_fatal("lc ct off!");}

#ifdef DEBUG
  error_msg_warning("det_rc() :: end w/%d of %d going right\n",nl,n);
#endif

  return(t_ptr);
}


/********************************sparse_matrix.c*******************************
Function: SMI_rsb_vertex()

Input : 
Output: 
Return: 
Description:  
*********************************sparse_matrix.c******************************/
int
SMI_rsb_vertex(int mat, int *color, int level, int *in_vertex, int nv)
{
  matrix_id_ptr m;
  int i, depth, rank, lda, *iptr, *order;
  stack_ADT s;
  queue_ADT in_q, out_q, tmp_q;
  struct tree_node *t_ptr, *t_lc, *t_rc;
  int *tmp_color;
  int nl, nr, ns, j, k;
  double *dptr, *rhs;
  double *fill_level, dtmp;
  double wt_subg_old=0.0;
  double wt_cuts_old=0.0;
  int *ia, *ja, n;
  double *vals;
  int *list, *vertex_order, nvu, nvo;
  int *vertex;
  double weight_old, weight_new;
  int  cuts_old, cuts_new;
  double *rhs2;
  int sq;
  int c, nc_start;
  int cost;
  double res, apl;
  int total_sep;


  /* check valid matrix handle */
  if ((mat<1)||(mat>MAX_MAT)||((m=active_matrix=matrix_list[--mat]) == NULL))
    {error_msg_fatal("SMI_rsb_vertex() :: matrix handle is NULL!!!");}

  /* problem size */
  lda = m->lda;
  start_nnz = nnz_adj_list(m->adj_list,lda);

#ifdef SAFE  
  /* check whether matrix is symmetric */
  if (!is_adj_list_symm(m->adj_list,lda))
    {error_msg_fatal("SMI_rsb_vertex() :: not symmetric!!!");}
#endif

  /* convert linked list into vector representation */
  compress_adj_list();

  /* queues/stack for building/processing tree */
  done_q = new_queue();  
  in_q   = new_queue();  
  out_q  = new_queue();  

  /* set tree root to original problem */
  /* later :: use -1 flag to indicate all lda rows */
  root = t_ptr = new_tree_node();
  t_ptr->num_act = lda;
  t_ptr->active = iptr = (int *) bss_malloc(INT_LEN*(lda+1));
  ivec_c_index(iptr,lda);
  *(iptr+lda) = -1;
  t_ptr->rhs = dptr = (double *) bss_malloc(2*lda*sizeof(double));
  Init_rhs_coord(dptr,iptr,lda);


  /* zero out vector and pre-label outflow bc vertices */
  /* 
  vertex = (int *) bss_malloc(INT_LEN*(nv+2));
  ivec_copy(vertex,in_vertex,nv+2);
  vertex_order = label_outflow(vertex,nv,&nvu,&nvo);
  */

  t_ptr->lb = 1;
  t_ptr->ub = nvu-nvo;

#ifdef DEBUG
  for (k=0; k<=nvu; k++)
    {printf("%d %d\n",k,vertex_order[k]);}
#endif

  /* queue root */
  enqueue(out_q,t_ptr);

  /* init color ==> elm to proc map) */
  ivec_set(color,-1,lda); 

  /* create binary tree and process all non-leaf nodes */
  weight_new = 1.0*start_nnz-lda;
  cuts_new = 0;
  cost = 0;
  max_depth = depth = MIN(log_2_floor(lda),level);

  printf("level 0 begin RSB w/nel=%d, nv=%d, nvu=%d, nvo=%d, dof=nvu-nvo=%d\n",
	 lda,nv,nvu,nvo,nvu-nvo);
  printf("level 0 begin RSB w/lanczos tolerance=%g, max repeats=%d, np max=%d\n",
	 res_limit,max_repeats,1<<max_depth);

  rsb = (struct rsb_info*) bss_malloc(sizeof(struct rsb_info));
  rsb->max_depth   = max_depth;
  rsb->depth       = depth;
  rsb->root        = root;
  rsb->nel         = lda;
  rsb->n           = nvu-nvo;
  rsb->done_q      = done_q;
  rsb->half_done_q = rsb->bad_q = NULL;

  total_sep=0;
  for(i=0; i<depth; i++)
    {
      weight_old = weight_new;
      weight_new = 0.0;
      cuts_old = cuts_new;
      cuts_new = 0;
      apl=0.0;
      for (rank=len_queue(out_q),j=0; j<rank; j++)
	{
	  /* subproblem parameters */
	  t_ptr = dequeue(out_q);
	  t_ptr->depth = i;
	  t_ptr->rank  = j;
	  n     = t_ptr->num_act;
	  list  = t_ptr->active;
	  rhs   = t_ptr->rhs;

	  /* extract corresponding submatrix */
	  weight_new += lap_sub_matrix(list, n, &ia, &ja, &vals);

	  /* dump_laplacian(n,ia,ja,vals); */

	  /* rsb on subgraph */
	  /* eventually place in sub routine */
	  tmp_color = sep_graph(list, n, ia, ja, vals, rhs, &nl);
	  k=0;
	  nc_start = nl;
	  res = res_limit;
	  if ((n>4)||(!tmp_color))
	    {
	      if (rsb_dim==_2D)
		{
		  sq = sqrt(1.0*n)+0.1;
		  c = 3*sq-2;
		  if (sq*sq - n)
		    {nr = max_repeats;}
		  else if (nl/2 > c)
		    {nr = 5*max_repeats;}
		  else
		    {nr = 2*max_repeats;}
		  res_limit = 1.0e-02;
		  for (k=0; ;k++)
		    {
		      res_limit = MAX(res_limit,1.0e-10);
		      hmt_perturb(rhs,n);
		      iptr = sep_graph(list, n, ia, ja, vals, rhs, &ns);
		      if (ns<nl)
			{
			  nl = ns;
			  if (tmp_color) 
			    {bss_free(tmp_color);}
			  tmp_color = iptr;
			}
		      else
			{
			  if (iptr)
			    {bss_free(iptr);}
			}

		      res_limit*=0.251189;

		      if (res_limit<1.0e-10)
			{res_limit = 1.0e-02;}

		      if ((k>nr)&&(tmp_color))
			{break;}

		      if (k>10*nr)
			{
			  error_msg_warning("utter sep_graph failure!\n");
			  break;
			}
		    }
		}
	      else
		{
		  sq = pow(1.0*n,1.0/3.0)+0.1;
		  c = 9*sq*sq-12*sq+4;
		  if (sq*sq*sq - n)
		    {nr = max_repeats;}
		  else if (nl/2 > c)
		    {nr = 5*max_repeats;}
		  else
		    {nr = 2*max_repeats;}
		  res_limit = 1.0e-02;
		  for (k=0; ; k++)
		    {
		      res_limit = MAX(res_limit,1.0e-10);
		      hmt_perturb(rhs,n);
		      iptr = sep_graph(list, n, ia, ja, vals, rhs, &ns);
		      if (ns<nl)
			{
			  nl = ns;
			  if (tmp_color)
			    {bss_free(tmp_color);}
			  tmp_color = iptr;
			}
		      else
			{
			  if (iptr)
			    {bss_free(iptr);}
			}

		      res_limit*=0.251189;

		      if (res_limit<1.0e-10)
			{res_limit = 1.0e-02;}

		      if ((k>nr)&&(tmp_color))
			{break;}

		      if (k>10*nr)
			{
			  error_msg_warning("utter sep_graph failure!\n");
			  break;
			}

		      hmt_perturb(rhs,n);
		    }
		}
	    }

	  apl+=1.0*k;
	  res_limit = res;

#ifdef DEBUG
	  printf("n=%d, #retrys=%d, sw=%d, fw=%d, c=%d,res=%g\n",n,k,
			    nc_start,nl,c,res_limit);
#endif

	  if (tmp_color==NULL)
	    {
	      error_msg_warning("SMI_rsb_vertex() :: tmp_color is NULL!\n");
	      break;
	    }
	  
	  if (nl==INT_MAX)
	    {
	      error_msg_fatal("SMI_rsb_vertex() :: nl eq. INT_MAX?\n");
	      break;
	    }

	  cuts_new += nl;

	  /* determine separator vertices */
	  nr = t_ptr->ub;
	  t_ptr->sep = det_sep(list,tmp_color,n,vertex,nv,vertex_order,nvu,nr,&ns);
	  t_ptr->num_sep   = ns;
	  total_sep += ns;

	  cost += ((t_ptr->ub-t_ptr->lb+1)*ns);

	  /* set up for rsb on left subgraph */
	  nr = t_ptr->lb;
	  t_ptr->lc = det_lc_vertex(list,tmp_color,rhs,n,vertex,nv,
				    vertex_order,nvu,nr,&nl);
	  enqueue(in_q,t_ptr->lc);

	  /* set up for rsb on right subgraph */
	  t_ptr->rc = det_rc_vertex(list,tmp_color,rhs,n,vertex,nv,
				    vertex_order,nvu,t_ptr->lb+nl,&nr);
	  enqueue(in_q,t_ptr->rc);

	  /* hold the entire tree as we build it */
	  enqueue(done_q,t_ptr);

	  /* might want to keep these for later use */
	  bss_free(tmp_color);
	  bss_free(ia);
	  bss_free(ja);
	  bss_free(vals);

	  nl = t_ptr->lc->num_act;	  
	  nr = t_ptr->rc->num_act;

	  /* missing an element? */
	  if (n!=(nl+nr))
	    {error_msg_fatal("SMI_rsb_vertex() :: (%d) ==> (%d,%d)\n",n,nl,nr);}
	}

      /* average number of retrys this level */
      apl/=1.0*j;      

      /* print out level info */
      printf("level %d done w/%d cuts, split %g==>(%g,%d) and apl=%.1f\n",
	     i,cuts_new,weight_old,weight_new,cuts_old,apl);

      /* if we're missing someone let us know */
      if (fabs(weight_old-(weight_new+1.0*cuts_old))>0.5)
	{
	  printf("again level %d done w/%d cuts and split %g ==> (%g,%d)\n",
		 i,cuts_new,weight_old,weight_new,cuts_old);
	}


      /* did we process the entire level? */
      if (!tmp_color)
	{
	  error_msg_warning("SMI_rsb_vertex() :: FAILURE!!!\n");
	  error_msg_warning("SMI_rsb_vertex() :: only made it to level %d\n",i);
	  error_msg_warning("SMI_rsb_vertex() :: Trying to recover ...\n");
	  fflush(stdout);

	  /* hold that bad apple */
	  rsb->bad_q = bad_q = new_queue();
	  enqueue(bad_q,t_ptr);

	  /* in_q now holds partially completed level */
	  rsb->depth = max_depth = i;
	  while (len_queue(out_q))
	    {
	      t_ptr = dequeue(out_q);
	      enqueue(in_q,t_ptr);
	    }
	  rsb->half_done_q = half_done_q = in_q;
	  break;
	}

      /* otherwise go on to next level */
      tmp_q = out_q;
      out_q = in_q;
      in_q  = tmp_q;
    }

  /* processed all levels */
  if (i==depth)
    {
      n = nvu-nvo;
      if (rsb_dim==_2D)
	{
	  printf("level %d done w/nnz = %d = %g * %d^(3/2)\n",
		 i,cost,cost/pow(1.0*n,1.5),n);
	}
      else
	{
	  printf("level %d done w/nnz = %d = %g * %d^(5/3)\n",
		 i,cost,cost/pow(1.0*n,5.0/3.0),n);
	}

      /* label leaves (stored in in_q) and queue on done_q */
      for (rank=len_queue(out_q),j=0; j<rank; j++)
	{
	  t_ptr = dequeue(out_q);
	  t_ptr->depth = i;
	  t_ptr->rank  = j;
	  n     = t_ptr->num_act;
	  list  = t_ptr->active;

#ifdef DEBUG
	  printf("i=%d, j=%d, lb=%d, ub=%d, n=%d, who=%d\n",i,j,
		 t_ptr->lb,t_ptr->ub,n,list[0]);
#endif
	  /* process */
	  nl = t_ptr->lb;
	  nr = t_ptr->ub;
	  if (nr<nl)
	    {
	      t_ptr->sep = NULL;
	      t_ptr->num_sep   = 0;
	    }
	  else
	    {
	      ns = 8;
	      t_ptr->sep = det_sep(list,NULL,n,vertex,nv,vertex_order,nvu,nr,&ns);
	      t_ptr->num_sep   = ns;
	      total_sep += ns;
	    }
	  enqueue(done_q,t_ptr);
	}
      free_queue(out_q);

      if (total_sep!=(nvu-nvo))
	{
	  printf("ts=%d, nvu=%d, nvo=%d ==> missing=%d\n",
		 total_sep,nvu,nvo,total_sep-(nvu-nvo));
	}
      else
	{
	  printf("ts=%d, nvu=%d, nvo=%d ==> none missing\n",
		 total_sep,nvu,nvo);
	}

      /* free tree */
      if (len_queue(in_q))
	{printf("SMI_rsb_vertex() :: in_q not empty!\n");}
      else
	{free_queue(in_q);}

      return(SUCCESS);
    }
  else
    {
      printf("level %d != %d ==> RSB Failure!\n",i,depth);
      free_queue(out_q);

      /* free tree */
      if (len_queue(in_q))
	{printf("SMI_rsb_vertex() :: in_q not empty!\n");}
      else
	{free_queue(in_q);}
      return(FAILURE);
    }
      
  /*
  bss_free(fill_level);
  bss_free(order);
  bss_free(vertex);
  bss_free(vertex_order);
  */
}



/********************************sparse_matrix.c*******************************
Function: SMI_rsb()

Input : 
Output: 
Return: 
Description:  
*********************************sparse_matrix.c******************************/
int
SMI_rsb(int mat, int *color, int level)
{
  matrix_id_ptr m;
  int i, depth, rank, lda, *iptr, *order;
  stack_ADT s;
  queue_ADT in_q, out_q, tmp_q;
  struct tree_node *t_ptr, *t_lc, *t_rc;
  int *tmp_color;
  int nl, nr, ns, j, k;
  double *dptr, *rhs;
  double *fill_level, dtmp;
  double wt_subg_old=0.0;
  double wt_cuts_old=0.0;
  int *ia, *ja, n;
  double *vals;
  int *list;
  double weight_old, weight_new;
  int  cuts_old, cuts_new;
  double *rhs2;
  int sq;
  int c, nc_start;
  int cost;
  double res, apl;
  int total_sep;

#ifdef DEBUG
  error_msg_warning("SMI_rsb() :: begin\n");
#endif

  /* check valid matrix handle */
  if ((mat<1)||(mat>MAX_MAT)||((m=active_matrix=matrix_list[--mat]) == NULL))
    {error_msg_fatal("SMI_rsb() :: matrix handle is NULL!!!");}

  /* problem size */
  lda = m->lda;
  start_nnz = nnz_adj_list(m->adj_list,lda);

#ifdef SAFE
  /* check whether matrix is symmetric */
  if (!is_adj_list_symm(m->adj_list,lda))
    {error_msg_fatal("SMI_rsb() :: not symmetric!!!");}
#endif

  /* convert linked list into vector representation */
  compress_adj_list();

  /* queues/stack for building/processing tree */
  done_q = new_queue();  
  in_q   = new_queue();  
  out_q  = new_queue();  

  /* set tree root to original problem */
  /* later :: use -1 flag to indicate all rows, 1...lda, active */
  root = t_ptr = new_tree_node();
  t_ptr->num_act = lda;
  t_ptr->active = iptr = (int *) bss_malloc(INT_LEN*(lda+1));
  ivec_c_index(iptr,lda);
  *(iptr+lda) = -1;
  /* why 2x on size? */
  t_ptr->rhs = dptr = (double *) bss_malloc(2*lda*sizeof(double));
  
  /* need seed vector ... use Cuthill-Mckee so build (i,j) pair list */
  /* should move this code into a subroutine ... later */
  n = lda;
  weight_new = lap_sub_matrix(t_ptr->active, n, &ia, &ja, &vals);
  i = ia[n] - n;
  if (i&1)
    {error_msg_fatal("SMI_rsb() :: matrix not symmetric :: %d!!\n",i);}
  if (i>0)
    {i>>=1;}
  iptr = order     = (int *) bss_malloc(INT_LEN*2*(i+1));
  list = tmp_color = (int *) bss_malloc(INT_LEN*n);
  *iptr++ = n;
  *iptr++ = k = i;
  for (i=j=0; i<n; i++)
    {
      while (j<ia[i+1])
	{
	  if (ja[j]>i)
	    {*iptr++ = i+1; *iptr++ = ja[j]+1; k--;}
	  j++;
	}
    }
  if (k)
    {error_msg_fatal("SMI_rsb() :: pair list failure!\n");}

  /* use min degree ordering */
  k=sm_bandwidth_reduction(order,tmp_color,0);
  if (k==0)
    {
      error_msg_warning("SMI_rsb() :: sm failure! Using random rhs init routine!\n");
      Init_rhs_random(dptr,lda);
    }
  else

    {
      if (k<0) {error_msg_warning("SMI_rsb() :: sm disconnected ... proceeding anyway!\n");}
      for (i=0;i<n;i++)
	{*dptr++ = *list++/(n+1.0);}
    }
  bss_free(ia);
  bss_free(ja);
  bss_free(vals);
  bss_free(order);
  bss_free(tmp_color);
  
  /* queue root */
  enqueue(out_q,t_ptr);

  /* init color ==> elm to proc map) */
  ivec_set(color,-1,lda); 

  /* create binary tree and process all non-leaf nodes */
  weight_new = 1.0*start_nnz-lda;
  cuts_new = 0;
  cost = 0;
  max_depth = depth = MIN(log_2_floor(lda),level);

  printf("level 0 begin RSB w/nel=%d\n",lda);
  printf("level 0 begin RSB w/lanczos tolerance=%g, max repeats=%d, np max=%d\n",
	 res_limit,max_repeats,1<<max_depth);
  fflush(stdout);
  

  rsb = (struct rsb_info*) bss_malloc(sizeof(struct rsb_info));
  rsb->max_depth   = max_depth;
  rsb->depth       = depth;
  rsb->root        = root;
  rsb->nel         = lda;
  rsb->done_q      = done_q;
  rsb->half_done_q = rsb->bad_q = NULL;

  total_sep=0;
  while(len_queue(out_q))
    {
      weight_old = weight_new;
      weight_new = 0.0;
      cuts_old = cuts_new;
      cuts_new = 0;
      apl=0.0;
      for (rank=len_queue(out_q),j=0; j<rank; j++)
	{
	  /* subproblem parameters */
	  t_ptr = dequeue(out_q);
	  t_ptr->depth = i;
	  t_ptr->rank  = j;
	  n     = t_ptr->num_act;
	  list  = t_ptr->active;
	  rhs   = t_ptr->rhs;

	  /* extract corresponding submatrix */
	  weight_new += lap_sub_matrix(list, n, &ia, &ja, &vals);

	  dump_laplacian(n,ia,ja,vals);

	  /* rsb on subgraph */
	  /* eventually place in sub routine */
	  tmp_color = sep_graph(list, n, ia, ja, vals, rhs, &nl);
	  k=0;
	  nc_start = nl;
	  res = res_limit;
	  if ((n>4)||(!tmp_color))
	    {
	      nr = 5;
	      res_limit = 1.0e-04;
	      for (k=0; ;k++)
		{
		  res_limit = MAX(res_limit,1.0e-10);
		  hmt_perturb(rhs,n);
		  iptr = sep_graph(list, n, ia, ja, vals, rhs, &ns);
		  if (ns<nl)
		    {
		      nl = ns;
		      if (tmp_color) 
			{bss_free(tmp_color);}
		      tmp_color = iptr;
		    }
		  else
		    {
		      if (iptr)
			{bss_free(iptr);}
		    }

		  res_limit*=0.251189;
		  if (res_limit<1.0e-10)
		    {res_limit = 1.0e-04;}

		  if ((k>nr)&&(tmp_color))
		    {break;}

		  if (k>10*nr)
		    {
		      error_msg_warning("utter sep_graph failure!\n");
		      break;
		    }
		}
	    }

	  apl+=1.0*k;
	  res_limit = res;

#ifdef DEBUG
	  printf("n=%d, #retrys=%d, sw=%d, fw=%d, c=%d,res=%g\n",n,k,
			    nc_start,nl,c,res_limit);
#endif

	  if (tmp_color==NULL)
	    {
	      printf("SMI_rsb() :: tmp_color is NULL!\n");
	      break;
	    }
	  
	  if (nl==INT_MAX)
	    {
	      error_msg_fatal("SMI_rsb() :: nl eq. INT_MAX?\n");
	      break;
	    }

	  cuts_new += nl;

	  /* determine separator vertices */
	  /*
	  nr = t_ptr->ub;
	  t_ptr->sep = det_sep(list,tmp_color,n,vertex,nv,vertex_order,nvu,nr,&ns);
	  t_ptr->num_sep   = ns;
	  total_sep += ns;

	  cost += ((t_ptr->ub-t_ptr->lb+1)*ns);
	  */

	  /* set up for rsb on left subgraph */
	  ivec_dump(tmp_color,n,n,n,"tmp_color before det_lc");	
	  t_ptr->lc = det_lc(list,tmp_color,rhs,n);
	  enqueue(in_q,t_ptr->lc);

	  /* set up for rsb on right subgraph */
	  t_ptr->rc = det_rc(list,tmp_color,rhs,n);
	  enqueue(in_q,t_ptr->rc);

	  /* hold the entire tree as we build it */
	  enqueue(done_q,t_ptr);

	  /* might want to keep these for later use */
	  bss_free(tmp_color);
	  bss_free(ia);
	  bss_free(ja);
	  bss_free(vals);

	  nl = t_ptr->lc->num_act;	  
	  nr = t_ptr->rc->num_act;

	  /* missing an element? */
	  if (n!=(nl+nr))
	    {error_msg_fatal("SMI_rsb() :: (%d) ==> (%d,%d)\n",n,nl,nr);}
	}

      /* average number of retrys this level */
      apl/=1.0*j;      

      /* print out level info */
      printf("level %d done w/%d cuts, split %g==>(%g,%d) and apl=%.1f\n",
	     i,cuts_new,weight_old,weight_new,cuts_old,apl);

      /* if we're missing someone let us know */
      if (fabs(weight_old-(weight_new+1.0*cuts_old))>0.5)
	{
	  printf("again level %d done w/%d cuts and split %g ==> (%g,%d)\n",
		 i,cuts_new,weight_old,weight_new,cuts_old);
	}


      /* did we process the entire level? */
      if (!tmp_color)
	{
	  printf("SMI_rsb() :: FAILURE!!!\n");
	  printf("SMI_rsb() :: only made it to level %d\n",i);
	  printf("SMI_rsb() :: Trying to recover ...\n");
	  fflush(stdout);

	  /* hold that bad apple */
	  rsb->bad_q = bad_q = new_queue();
	  enqueue(bad_q,t_ptr);

	  /* in_q now holds partially completed level */
	  rsb->depth = max_depth = i;
	  while (len_queue(out_q))
	    {
	      t_ptr = dequeue(out_q);
	      enqueue(in_q,t_ptr);
	    }
	  rsb->half_done_q = half_done_q = in_q;
	  break;
	}

      /* otherwise go on to next level */
      tmp_q = out_q;
      out_q = in_q;
      in_q  = tmp_q;
    }

  /* processed all levels */
  if (i==depth)
    {
      for (rank=len_queue(out_q),j=0; j<rank; j++)
	{
	  t_ptr = dequeue(out_q);
	  t_ptr->depth = i;
	  t_ptr->rank  = j;
	  n     = t_ptr->num_act;
	  list  = t_ptr->active;

#ifdef DEBUG
	  printf("i=%d, j=%d, lb=%d, ub=%d, n=%d, who=%d\n",i,j,
		 t_ptr->lb,t_ptr->ub,n,list[0]);
#endif
	  /* process */
	  t_ptr->sep = NULL;
	  t_ptr->num_sep   = 0;

	  enqueue(done_q,t_ptr);
	}
      free_queue(out_q);

      /* free tree */
      if (len_queue(in_q))
	{error_msg_warning("SMI_rsb() :: in_q not empty!\n");}
      else
	{free_queue(in_q);}

      return(SUCCESS);
    }
  else
    {
      error_msg_warning("level %d != %d ==> RSB Failure!\n",i,depth);
      free_queue(out_q);

      /* free tree */
      if (len_queue(in_q))
	{error_msg_warning("SMI_rsb() :: in_q not empty!\n");}
      else
	{free_queue(in_q);}
      return(FAILURE);
    }
      
  /*
  bss_free(fill_level);
  bss_free(order);
  bss_free(vertex);
  bss_free(vertex_order);
  */
}



/********************************sparse_matrix.c*******************************
Function: SMI_()

Input : 
Output: 
Return: 
Description:  
*********************************sparse_matrix.c******************************/
void
rsb_peep_ ()
{
  int i,j,k;
  int nl, nq,na,off;
  int *list;
  queue_ADT tmp_q, hold_q;
  struct tree_node *t_ptr;
  int color_jump;


  if (bad_q)
    {
      error_msg_warning("rsb_peep_() :: RSB Failed here!\n");

      if (len_queue(bad_q))
	{
	  t_ptr = (struct tree_node *) dequeue(bad_q);
	  enqueue(bad_q,t_ptr);
	  
	  na = t_ptr->num_act;
	  list = t_ptr->active;
	  error_msg_warning("RSB choked on the following %d elements :: ",na);
	  
	  for (i=0; i<na; i++)
	    {
	      j = list[i]+1;
	      k=2;
	      /* hmt_drawel_color_(&j,&k); */
	      printf("%d ",list[i]+1);
	    }
	  return;
	}
    }
  else
    {

      error_msg_warning("rsb_peep() :: RSB didn't fail!$");
      return;
    }    
}


/********************************sparse_matrix.c*******************************
Function: SMI_()

Input : 
Output: 
Return: 
Description:  
*********************************sparse_matrix.c******************************/
void
SMI_extract_proc_map(int *color, int level)
{
  int i,j,k;
  int nl, nq,na,off;
  int *list;
  queue_ADT tmp_q, hold_q;
  struct tree_node *t_ptr;
  int color_jump;
  int num_colors;


  /* num_colors = x_num_colors(); */
  num_colors = 1;

  hold_q = new_queue();

  level = MIN(level,max_depth);

  nl = 1;
  i = level;
  while (i--)
    {nl<<=1;}

  if (nl>num_colors)
    {
      error_msg_warning("don't have that many colors!\n");
      color_jump = 1;
    }
  else
    {color_jump = num_colors/nl;}

  if (level==0)
    {off=0;}
  else
    {off = (1<<level)-1;}

  nq = len_queue(done_q);
#ifdef DEBUG
  printf("level=%d, nl=%d, off=%d, nq=%d\n",level,nl,off,nq);
#endif
  if ((off+nl) > nq)
    {
      if (bad_q)
	{
	  error_msg_warning("RSB Failed here!\n");

	  if (len_queue(bad_q))
	    {
	      t_ptr = (struct tree_node *) dequeue(bad_q);
	      enqueue(bad_q,t_ptr);
	      
	      na = t_ptr->num_act;
	      list = t_ptr->active;
	      error_msg_warning("RSB choked on the following %d elements :: ",na);
	      ivec_set(color,-1,5000); 
	      for (j=0; j<na; j++)
		{
		  printf("%d ",list[j]+1);
		  color[list[j]] = 2;
		}
	      return;
	    }
	}

#ifdef DEBUG
      error_msg_warning("SMI_extract_proc() :: not that many in queue!\n");
#endif

      error_msg_warning("Can't show that many procs!\n");
      return;
    }    
  
  for (i=0; i<off; i++)
    {
      t_ptr = (struct tree_node *) dequeue(done_q);
      enqueue(hold_q,t_ptr);
    }

  nl += off;
  for (k=1,i=off; i<nl; i++,k++)
    {
      t_ptr = (struct tree_node *) dequeue(done_q);
      enqueue(hold_q,t_ptr);
      
      na = t_ptr->num_act;
      list = t_ptr->active;
#ifdef DEBUG      
      printf("p%4d has %4d elms :: ",k,na);
#endif
      for (j=0; j<na; j++)
	{
#ifdef DEBUG      
	  printf("%d ",list[j]+1);
#endif
	  color[list[j]] = (k*color_jump)%num_colors+16;
	}
#ifdef DEBUG      
      printf("\n");
#endif
    }
#ifdef DEBUG      
  printf("\n");
#endif

  for (i=nl; i<nq; i++)
    {
      t_ptr = (struct tree_node *) dequeue(done_q);
      enqueue(hold_q,t_ptr);
    }

  free_queue(done_q);
  rsb->done_q = done_q = hold_q;
}



/********************************sparse_matrix.c*******************************
Function: SMI_()

Input : 
Output: 
Return: 
Description:  
*********************************sparse_matrix.c******************************/
void
SMI_dump_xxt(char *file_prefix, int length)
{
  int i,j,k,ct;
  int nl, ns, nv, nc, n, nel, depth;
  int tns=0, tnl=0, total, rank;
  int *list, *sep, *vertex, *vertex_map, *remapped_vertex;
  int *proc, *elms;
  int *iptr;
  queue_ADT q, hold_q;
  struct tree_node *t_ptr;
  FILE *ofp_map, *ofp_sep;
  char name1[C_LINE+1];
  char name2[C_LINE+1];


  /* did rsb go alll the way out? */
  if (!rsb)
    {
      error_msg_warning("You have to \"DO IT\" first!\n");
      return;
    }

  if (rsb->max_depth!=rsb->depth)
    {error_msg_warning("can't dump because tree not deep enough!\n"); return;}

#ifdef OLD
  if ((ofp_map=fopen("xxt_map.rea","w"))==NULL)
    {error_msg_fatal("can't open xxt_map for output!\n"); return;}

  if ((ofp_sep=fopen("xxt_sep.rea","w"))==NULL)
    {error_msg_fatal("can't open xxt_sep for output!\n"); return;}
#endif

  if (length>C_LINE)
    {error_msg_fatal("SMI_dump_xxt():: rea prefix too long!\n"); return;}

  for (i=0;i<length;i++) 
      {
       name1[i]=file_prefix[i];
      }
  name1[i]='\0';
  strcat(name1,".map");

  for (i=0;i<length;i++) 
      {
       name2[i]=file_prefix[i];
      }
  name2[i]='\0';
  strcat(name2,".sep");

  if ((ofp_map=fopen(name1,"w"))==NULL)
    {error_msg_fatal("SMI_dump_xxt():: map file %s not found!\n",name1);return;}

  if ((ofp_sep=fopen(name2,"w"))==NULL)
    {error_msg_fatal("SMI_dump_xxt():: sep file %s not found!\n",name2);return;}

  /* problem parameters */
  depth  = rsb->depth;
  nel = rsb->nel;
  nv  = rsb->nv;
  n   = rsb->nvu-rsb->nvo;
  nc  = (rsb_dim==_2D) ? (4) : (8);
  if ((nc*nel)!=nv)
    {error_msg_fatal("nv != nel*nc\n");}
  
  /* these become file headers */
  fprintf(ofp_map,"%-8d %-8d %-8d %-8d %-8d %-8d %-8d\n",nel,n,depth,
	  1<<depth,nv,rsb->nvu,rsb->nvo);
  fprintf(ofp_sep,"%-8d %-8d %-8d %-8d %-8d %-8d %-8d\n",nel,n,depth,
	  1<<depth,nv,rsb->nvu,rsb->nvo);

  /* remap vertices and shift to hypercube frame of reference */
  vertex     = rsb->vertex;
  vertex_map = rsb->vertex_map;
  remapped_vertex = iptr = (int *) bss_malloc(nv*INT_LEN);
  nc>>=1;
  for (k=i=0; i<nel; i++)
    {
      for (j=0; j<nc; j++)
	{
	  if (j%2)
	    {
	      *iptr++ = vertex_map[vertex[k+1]];
	      *iptr++ = vertex_map[vertex[k]];
	    }
	  else
	    {
	      *iptr++ = vertex_map[vertex[k]];
	      *iptr++ = vertex_map[vertex[k+1]];
	    }
	  k+=2;
	}
    }
  nc<<=1;
  if (k!=nv)
    {error_msg_fatal("k != nv\n");}

  proc  = (int *) bss_malloc(nel*INT_LEN);
  elms  = (int *) bss_malloc(nel*INT_LEN);


  q      = rsb->done_q;
  hold_q = new_queue();
  total = len_queue(q);

  if (total!=((1<<(depth+1))-1))
    {error_msg_fatal("lenq=%d, depth=%d\n",total,depth);}
  
  tnl = tns = ct = 0;
  for (i=0; i<depth; i++)
    {
      rank = 1<<i;
      for (j=0; j<rank; j++)
	{
	  if (ct++==total)
	    {error_msg_fatal("again hmt can't count!\n");}

	  enqueue(hold_q,t_ptr=dequeue(q));
	  
	  sep = t_ptr->sep;
	  tns += ns = t_ptr->num_sep;
	  fprintf(ofp_sep,"%-4d ",ns);
	  for (k=0; k<ns; k++)
	    {fprintf(ofp_sep,"%-6d ",*sep++);}
	  fprintf(ofp_sep,"\n");
	}
    }

  for ( ; i<=depth; i++)
    {
      rank = 1<<i;
      
      for (j=0; j<rank; j++)
	{
	  if (ct++==total)
	    {error_msg_fatal(" again again hmt can't count!\n");}

	  enqueue(hold_q,t_ptr=dequeue(q));

	  sep = t_ptr->sep;
	  tns += ns = t_ptr->num_sep;
	  fprintf(ofp_sep,"%-4d ",ns);
	  for (k=0; k<ns; k++)
	    {fprintf(ofp_sep,"%-6d ",*sep++);}
	  fprintf(ofp_sep,"\n");

	  list = t_ptr->active;
	  nl = t_ptr->num_act;
	  if ((nl!=1)&&(nl!=2))
	    {error_msg_fatal("1 or 2 elms per proc not %d!\n",nl);}

	  for (k=0; k<nl; k++)
	    {proc[tnl] = j; elms[tnl] = (*list++ + 1); tnl++;}
	}
    }
  if (tnl != nel)
    {error_msg_fatal("nel=%d but number found=%d!\n",nel,tnl);}

  if (tns != n)
    {error_msg_fatal("n=%d but number found=%d!\n",n,tns);}

  /* go from proc elm to elm proc map */
  ivec_sort_companion(elms,proc,nel);
  iptr = remapped_vertex;
  for (i=0; i<nel; i++)
    {
      fprintf(ofp_map,"%-4d ",proc[i]);
      for (k=0; k<nc; k++)
	{fprintf(ofp_map,"%-6d ", *iptr++);}
      fprintf(ofp_map,"\n");	  
    }

  fclose(ofp_map);
  fclose(ofp_sep);

  bss_free(proc);
  bss_free(elms);
  bss_free(remapped_vertex);
	  
  free_queue(q);
  rsb->done_q = done_q = hold_q;
}



/********************************sparse_matrix.c*******************************
Function: SMI_()

Input : 
Output: 
Return: 
Description:  
*********************************sparse_matrix.c******************************/
void
clear_rsb_peep_ ()
{
  struct tree_node *t_ptr;
  queue_ADT q;


  if (q=rsb->bad_q)
    {
      while(len_queue(q))
	{
	  t_ptr = dequeue(q);
	  if (t_ptr->num_act)  {bss_free(t_ptr->active);}
	  if (t_ptr->sep)      {bss_free(t_ptr->sep);}
	  if (t_ptr->rhs)      {bss_free(t_ptr->rhs);}
	  bss_free(t_ptr);
	}
      free_queue(q);
      rsb->bad_q=bad_q=NULL;
    }
}



/********************************sparse_matrix.c*******************************
Function: SMI_()

Input : 
Output: 
Return: 
Description:  
*********************************sparse_matrix.c******************************/
SMI_rsb_clean()
{
  struct tree_node *t_ptr;
  queue_ADT q;

  root=NULL;

  if (rsb)
    {
      if (rsb->vertex)     {bss_free(rsb->vertex);}
      if (rsb->vertex_map) {bss_free(rsb->vertex_map);}

      if (q=rsb->done_q)
	{
	  while(len_queue(q))
	    {
	      t_ptr = dequeue(q);
	      if (t_ptr->num_act)  {bss_free(t_ptr->active);}
	      if (t_ptr->sep)      {bss_free(t_ptr->sep);}
	      if (t_ptr->rhs)      {bss_free(t_ptr->rhs);}
	      bss_free(t_ptr);
	    }
	  free_queue(q);
	  done_q = NULL;
	}

      if (q=rsb->half_done_q)
	{
	  while(len_queue(q))
	    {
	      t_ptr = dequeue(q);
	      if (t_ptr->num_act)  {bss_free(t_ptr->active);}
	      if (t_ptr->sep)      {bss_free(t_ptr->sep);}
	      if (t_ptr->rhs)      {bss_free(t_ptr->rhs);}
	      bss_free(t_ptr);
	    }
	  free_queue(q);
	  half_done_q=NULL;
	}  

      if (q=rsb->bad_q)
	{
	  while(len_queue(q))
	    {
	      t_ptr = dequeue(q);
	      if (t_ptr->num_act)  {bss_free(t_ptr->active);}
	      if (t_ptr->sep)      {bss_free(t_ptr->sep);}
	      if (t_ptr->rhs)      {bss_free(t_ptr->rhs);}
	      bss_free(t_ptr);
	    }
	  free_queue(q);
	  bad_q=NULL;
	}

      bss_free(rsb);
      rsb=NULL;
    }
}



/********************************sparse_matrix.c*******************************
Function: hmt_ortho()

Input : 
Output: 
Return: 
Description:  
*********************************sparse_matrix.c******************************/
void 
hmt_ortho(double *r, int n)
{
  int i;
  double sum = 0.0;


  for (i=0; i<n; i++)
    {sum+=r[i];}

  sum/=n;

  for (i=0; i<n; i++)
    {r[i] -= sum;}
}



/********************************sparse_matrix.c*******************************
Function: normalize()

Input : 
Output: 
Return: 
Description:  
*********************************sparse_matrix.c******************************/
void 
normalize(double *r, int n)
{
  double norm;

  
  norm = 1.0/sqrt(hmt_ddot(r,r,n));
  while (n--)
    {*r++ *= norm;}
}



/********************************sparse_matrix.c*******************************
Function: Init_guess()

Input : 
Output: 
Return: 
Description:  seed x_0 for CG. Take it to be zero.
*********************************sparse_matrix.c******************************/
double *Init_guess(int n)
{
  double *rhs;


  rhs = (double *) bss_malloc(n*sizeof(double));
  dvec_zero(rhs,n);
  return(rhs);
}



/********************************sparse_matrix.c*******************************
Function: hmt_ddot() 

Input : 
Output: 
Return: 
Description:  vector vector dot product: <x1,x2> = x1^t.x2
*********************************sparse_matrix.c******************************/
double hmt_ddot(register double *x1, register double *x2, register int n)
{
  register double tmp=0.0;


  while (n--) {tmp+= (*x1++ * *x2++);}
  return(tmp);
}



/********************************sparse_matrix.c*******************************
Function: hmt_daxpy()

Input : 
Output: 
Return: 
Description:   adds two vectors: x1+x2
*********************************sparse_matrix.c******************************/
void
hmt_daxpy(register double *x1, register double scale, register double *x2, 
	  register int n)
{
  while (n--) {*x1++ += (scale * *x2++);}
}



/********************************sparse_matrix.c*******************************
Function: rev_daxpy()

Input : 
Output: 
Return: 
Description:   adds two vectors: x1+x2
*********************************sparse_matrix.c******************************/
void
rev_daxpy(register double *x1, register double scale, register double *x2, 
	  register int n)
{
  while (n--) {*x2 *= scale; *x2++ += *x1++;}
}



/********************************sparse_matrix.c*******************************
Function: matrix_mult()

Input : 
Output: 
Return: 
Description:   multiplies a vector by a real constant: alph*ax2_i.
*********************************sparse_matrix.c******************************/
void
matrix_mult(register double *Ap_k, register double *p_k, register int n, 
	    register int *ia, register int *ja, register double *vals)
{
  register int i;

  for (i=0, ia++; *ia>0; Ap_k++, ia++)
    {
      *Ap_k = 0.0;
      while (i<*ia)
	{*Ap_k += (*vals * p_k[*ja]); vals++; ja++; i++;}
    }
}



/********************************sparse_matrix.c*******************************
Function: matrix_mult()

Input : 
Output: 
Return: 
Description:   multiplies a vector by a real constant: alph*ax2_i.
*********************************sparse_matrix.c******************************/
void
skew_matrix_mult(register double *Ap_k, register double *p_k, register int n, 
		 register int *ia, register int *ja, register double *vals, 
		 double lambda)
{
  register int i, j;

  j=0;
  for (i=0, ia++; *ia>0; Ap_k++, ia++, j++)
    {
      *Ap_k = 0.0;
      while (i<*ia)
	{
	  if (*ja==j)
	    {*Ap_k += ((*vals-lambda) * p_k[*ja]); vals++; ja++; i++;}
	  else
	    {*Ap_k += (*vals * p_k[*ja]); vals++; ja++; i++;}
	}
    }
}



void
break_symm_lap(int n, int *ia, int *ja, double *vals)
{
  register int i, j=0, k=0;
  double input_f;
  double f, rm;
  char buf[81];

  printf("Adj Diag: ");
  if (fgets(buf, sizeof(buf), stdin) == NULL)
    {error_msg_fatal("gets() failure in break_symm_lap()\n");}
  input_f = atof(buf);

  srand(1024);

  /*
  f = 0.001 + (1.0*n)/500;
  f = MIN(f,0.5);
  */
  f = (double) input_f;
  k=n%2;
  rm=0.0;
  for (i=0, ia++; *ia>0; ia++, j++)
    {
      while (i<*ia)
	{
	  if (*ja==j)
	    {
	      *vals += f * ((double) rand()/RAND_MAX);
	    }
	  ja++;
	  vals++;
	  i++;
	}
    }
  if (fabs(rm)>1.0e-14)
    {error_msg_warning("bad break symm!\n");}
}


void
fdump_laplacian(int n, int *ia, int *ja, double *vals)
{
  int i, j, ct=0;
  FILE *ofp;


  ia++;
  ofp = fopen("A","w");
  fprintf(ofp,"A = \n");
  fprintf(ofp,"{ \n");
  for (i=0; i<n; i++)
    {
      fprintf(ofp,"{ ");
      if (*ja==0)
	{fprintf(ofp,"%.1f",*vals); ja++; vals++; ct++;}
      else
	{fprintf(ofp,"%.1f",0.0);}
      for (j=1; j<n; j++)
	{
	  if ((*ja==j)&&(ct<*ia))
	    {fprintf(ofp,", %.1f",*vals); ja++; vals++; ct++;}
	  else
	    {fprintf(ofp,", %.1f",0.0);}
	}
      if (i!=n-1)
	{fprintf(ofp,"},\n");}
      else
	{fprintf(ofp,"}\n");}
      ia++;
    }
  fprintf(ofp,"}\n");
  fclose(ofp);
}




void
dump_laplacian(int n, int *ia, int *ja, double *vals)
{
  int i, j, ct=0;

  ia++;
  printf("{ \n");
  for (i=0; i<n; i++)
    {
      printf("{ ");
      if (*ja==0)
	{printf("%.1f",*vals); ja++; vals++; ct++;}
      else
	{printf("%.1f",0.0);}
      for (j=1; j<n; j++)
	{
	  if ((*ja==j)&&(ct<*ia))
	    {printf(", %.1f",*vals); ja++; vals++; ct++;}
	  else
	    {printf(", %.1f",0.0);}
	}
      if (i!=n-1)
	{printf("},\n");}
      else
	{printf("}\n");}
      ia++;
    }
  printf("}\n");
  fflush(stdout);
}

void
full_dump_laplacian(int n, int *ia, int *ja, double *vals)
{
  int i, j, ct=0;

  ia++;
  printf("{ \n");
  for (i=0; i<n; i++)
    {
      printf("{ ");
      if (*ja==0)
	{printf("%.3f",*vals); ja++; vals++; ct++;}
      else
	{printf("%.1f",0.0);}
      for (j=1; j<n; j++)
	{
	  if ((*ja==j)&&(ct<*ia))
	    {printf(", %.3f",*vals); ja++; vals++; ct++;}
	  else
	    {printf(", %.1f",0.0);}
	}
      if (i!=n-1)
	{printf("},\n");}
      else
	{printf("}\n");}
      ia++;
    }
  printf("}\n");
  fflush(stdout);
}




void
diag_dump_laplacian(int n, int *ia, int *ja, double *vals)
{
  int i, j, ct=0;

  ia++;
  printf("{ \n");
  for (i=0; i<n; i++)
    {
      printf("{ ");
      if (*ja==0)
	{
	  if (*ja==i)
	    {printf("%.16f",*vals); ja++; vals++; ct++;}
	  else
	    {printf("%.1f",*vals); ja++; vals++; ct++;}
	}
      else
	{printf("%.1f",0.0);}
      for (j=1; j<n; j++)
	{
	  if ((*ja==j)&&(ct<*ia))
	    {
	      if (*ja==i)
		{printf(", %.16f",*vals); ja++; vals++; ct++;}
	      else
		{printf(", %.1f",*vals); ja++; vals++; ct++;}
	    }
	  else
	    {printf(", %.1f",0.0);}
	}
      if (i!=n-1)
	{printf("},\n");}
      else
	{printf("}\n");}
      ia++;
    }
  printf("}\n");
  fflush(stdout);
}


/********************************sparse_matrix.c*******************************
Function: lanczos()

Input : 
Output: 
Return: 
Description:  
*********************************sparse_matrix.c******************************/
double *
lanczos(double *rhs, int n, int *ia, int *ja, double *vals)
{
  int i,j,k=0;
  double *r_k, *p_k, *x_k, *Ap_k, *r_k_old, *r_k_new;
  double alpha=0.0, beta=0.0, rr=1.0, rr_limit=0.0, new_rr=0.0;
  double alpha_old=0.0, beta_old=0.0;
  int iter_limit;
  queue_ADT ln_q;
  struct lanczos_node *ln_ptr;
  double  *dv, *ev, *dv_w, *ev_w, **Q;
  double *temp;
  int p, b;
  double *new_rhs;

  /* queue to hold lanczos vectors */

  p = perm_calls()-perm_frees();
  b = bss_calls()-bss_frees();

  ln_q = new_queue();

  /* make sure rhs is orthogonal to (1,...,1)^T and has unit 2-norm */
  hmt_ortho(rhs,n);
  normalize(rhs,n);


  /* grab space and init */
  /*
  r_k     = (double *) bss_malloc((n+1)*sizeof(double));
  r_k_old = (double *) bss_malloc((n+1)*sizeof(double));
  r_k_new = (double *) bss_malloc((n+1)*sizeof(double));
  p_k     = (double *) bss_malloc((n+1)*sizeof(double));
  x_k     = (double *) bss_malloc((n+1)*sizeof(double));
  Ap_k    = (double *) bss_malloc((n+1)*sizeof(double));
  */
  new_rhs = (double *) bss_malloc((2*n+1)*sizeof(double));
  r_k     = (double *) bss_malloc((2*n+1)*sizeof(double));
  r_k_old = (double *) bss_malloc((2*n+1)*sizeof(double));
  r_k_new = (double *) bss_malloc((2*n+1)*sizeof(double));
  p_k     = (double *) bss_malloc((2*n+1)*sizeof(double));
  x_k     = (double *) bss_malloc((2*n+1)*sizeof(double));
  Ap_k    = (double *) bss_malloc((2*n+1)*sizeof(double));
  dvec_zero(r_k,2*n+1);
  dvec_zero(r_k_old,2*n+1);
  dvec_zero(r_k_new,2*n+1);
  dvec_zero(p_k,2*n+1);
  dvec_zero(x_k,2*n+1);
  dvec_zero(Ap_k,2*n+1);


  dvec_copy(r_k,rhs,n);
  dvec_copy(r_k_old,rhs,n);
  dvec_copy(p_k,rhs,n);
  dvec_zero(r_k_new,n);
  dvec_zero(x_k,n);
  dvec_zero(Ap_k,n);

  /* set lanczos tolerences */
  /* rr = fddot(n,r_k,1,r_k,1); */
  rr = hmt_ddot(r_k,r_k,n);
  rr_limit = res_limit;
  iter_limit = MIN(n,CG_LIMIT);


#ifdef DEBUG
  printf("\n");
  printf("rr %.16e\n",rr);
  printf("rr_limit %.16e\n",rr_limit);
  printf("iter_limit %d\n",iter_limit);
  printf("orthog_limit %.16e\n",RES_ORTHOG);
  printf("\n");
#endif

#ifdef MYMALLOC  
   perm_init(); 
#endif

  Q = (double **) perm_malloc((n+1)*sizeof(double *));

#ifdef DEBUG
  printf("rr[%4d]=%.16e\n",0,rr);
#endif

  for (;;)
    {
      k++;
      ln_ptr = (struct lanczos_node *)bss_malloc(sizeof(struct lanczos_node));
      ln_ptr->v = (double *) perm_malloc((n+1)*sizeof(double));
      dvec_copy(ln_ptr->v,r_k,n);
      normalize(ln_ptr->v,n);
      dvec_copy(r_k_old,ln_ptr->v,n);
      
      matrix_mult(Ap_k,p_k,n,ia,ja,vals);
      /* alpha = rr/fddot(n,Ap_k,1,p_k,1); */
      alpha = rr/hmt_ddot(Ap_k,p_k,n);
      /* fdaxpy(n,alpha,p_k,1,x_k,1); */
      hmt_daxpy(x_k,alpha,p_k,n);
      /* fdaxpy(n,-alpha,Ap_k,1,r_k,1); */
      hmt_daxpy(r_k,-alpha,Ap_k,n);
      hmt_ortho(r_k,n);

      dvec_copy(r_k_new,r_k,n);
      normalize(r_k_new,n);

      /* new_rr = fddot(n,r_k,1,r_k,1); */
      new_rr = hmt_ddot(r_k,r_k,n);
      beta = new_rr/rr;
      rev_daxpy(r_k,beta,p_k,n);
      /* hmt_ortho(p_k,n); */

      if (k==1)
	{
	  ln_ptr->d = 1.0/alpha; 
	  ln_ptr->e = 0.0;
	}
      else
	{
	  ln_ptr->d = beta_old/alpha_old + 1.0/alpha; 
	  ln_ptr->e = -sqrt(beta_old)/alpha_old;
	}
      enqueue(ln_q,ln_ptr);

      rr = new_rr;
      alpha_old = alpha;
      beta_old = beta;
    
#ifdef DEBUG      
      printf("2-norm r_k[%d]=%.16e\n",k,rr);
#endif


      if (rr<rr_limit)
	{
#ifdef DEBUG      
	  printf("rr[%4d]=%.16e\n",k,rr);
	  printf("rr_limit reached\n");
#endif
	  break;
	}

      if (rr>100.0)
	{
#ifdef DEBUG      
	  printf("rr[%4d]=%.16e\n",k,rr);
#endif
	  printf("rr_limit exceeded!!!\n");

	  bss_free(new_rhs);
	  bss_free(r_k);
	  bss_free(r_k_old);
	  bss_free(r_k_new);
	  bss_free(p_k);
	  bss_free(x_k);
	  bss_free(Ap_k);
	  perm_free(Q);

	  for (i=1; i<=k; i++)
	    {ln_ptr = dequeue(ln_q); bss_free(ln_ptr);}
	  
	  if (len_queue(ln_q))
	    {printf("lanczos queue not exhasuted!");}
	  else
	    {free_queue(ln_q);}
	  return(NULL);
	}

      if (k>=iter_limit)
	{
#ifdef DEBUG      
	  printf("rr[%4d]=%.16e\n",k,rr);
#endif
	  printf("iter_limit exceeded\n");
	  bss_free(new_rhs);
	  bss_free(r_k);
	  bss_free(r_k_old);
	  bss_free(r_k_new);
	  bss_free(p_k);
	  bss_free(x_k);
	  bss_free(Ap_k);
	  perm_free(Q);

	  for (i=1; i<=k; i++)
	    {ln_ptr = dequeue(ln_q); bss_free(ln_ptr);}
	  
	  if (len_queue(ln_q))
	    {printf("lanczos queue not exhasuted!");}
	  else
	    {free_queue(ln_q);}

	  return(NULL);
	}

      /* if (fabs(fddot(n,r_k_new,1,r_k_old,1)) >= RES_ORTHOG) */
      if (fabs(hmt_ddot(r_k_new,r_k_old,n)) >= RES_ORTHOG)
	{
#ifdef DEBUG      
	  printf("2-norm r_k[%4d]=%.16e\n",k,rr);
	  /* printf("(r[%4d],r[%4d])=%.16e \n",k,k-1,
	     fabs(fddot(n,r_k_new,1,r_k_old,1))); */
	  printf("(r[%4d],r[%4d])=%.16e \n",k,k-1,
		 fabs(hmt_ddot(r_k_new,r_k_old,n)));
#endif
	  printf("lost orthogonality of lanczos vectors\n");
	  bss_free(new_rhs);
	  bss_free(r_k);
	  bss_free(r_k_old);
	  bss_free(r_k_new);
	  bss_free(p_k);
	  bss_free(x_k);
	  bss_free(Ap_k);
	  perm_free(Q);

	  for (i=1; i<=k; i++)
	    {ln_ptr = dequeue(ln_q); bss_free(ln_ptr);}
	  
	  if (len_queue(ln_q))
	    {printf("lanczos queue not exhasuted!");}
	  else
	    {free_queue(ln_q);}
	  return(NULL);
	}
    }
#ifdef DEBUG
  if (k<2)
    {printf("Only one lanczos iteration?\n");}
#endif
  
  /* check cg convergence */
  matrix_mult(Ap_k,x_k,n,ia,ja,vals);
  dvec_copy(r_k,rhs,n);
  alpha=0.0;
  for (i=0; i<n; i++)
    {alpha = MAX(fabs(Ap_k[i]-r_k[i]),alpha);}

#ifdef DEBUG
  printf("inf-norm Ax_%3.d-b = %.16e\n",k,alpha);
#endif

  dv   = r_k_old;
  ev   = r_k_new;
  dv_w = p_k;
  ev_w = x_k;

  dvec_zero(dv,k+1);
  dvec_zero(dv_w,k+1);
  dvec_zero(ev,k+1);
  dvec_zero(ev_w,k+1);

  /*
  dvec_zero(dv,n+1);
  dvec_zero(dv_w,n+1);
  dvec_zero(ev,n+1);
  dvec_zero(ev_w,n+1);
  */

  for (i=1; i<=k; i++)
    {
      ln_ptr = dequeue(ln_q);
      Q[i] = ln_ptr->v;
      
      dv[i] = dv_w[i] = ln_ptr->d;
      if (i>1)
	{ev[i] =  ev_w[i] = ln_ptr->e;}
      bss_free(ln_ptr);
    }

  if (len_queue(ln_q))
    {printf("lanczos queue not exhasuted!");}
  else
    {free_queue(ln_q);}

  tqli(dv,ev,k,NULL);
  dvec_sort(dv,k+1);

#ifdef DEBUG
  for (i=0;i<=k;i++)
    {printf("ev[%d]=%.16e\n",i,dv[i]);}
#endif

  beta = dv[1];
#ifdef DEBUG      
  printf("min eval = %.16e\n",dv[1]);
  if (k>1)
    {printf("sec eval = %.16e\n",dv[2]);}
#endif    

#ifdef DEBUG
  printf("max eval = %.16e\n",dv[k]);
  printf("cond num = %.16e\n",dv[k]/dv[1]);
  printf("begin fiedler vector calc\n");
#endif

  /*  tri_calc_evec(dv_w,ev_w,k,beta,r_k); */
  r_k[0]=0.0;
  dvec_set(r_k+1,1.0,k);
  normalize(r_k+1,k);
  j=0;
  temp=trid_mult(dv_w,ev_w,k,r_k);
  alpha=0.0;
  for (i=1; i<=k; i++)
    {alpha = MAX(alpha,fabs(beta*r_k[i]-temp[i]));}
  bss_free(temp);
  /* printf("trid[%d] inf-norm Aev-lev = %.16e\n",j,alpha); */
  rr=0.0;
  for (;;)
    {    
      j++;
      if (trid_slv(dv_w,ev_w,k,beta-rr,r_k))
	{
	  temp = trid_mult(dv_w,ev_w,k,r_k);
	  alpha=0.0;
	  for (i=1; i<=k; i++)
	    {alpha = MAX(alpha,fabs(beta*r_k[i]-temp[i]));}

	  if ((alpha<1.0e-15)||(j))
	    {
#ifdef DEBUG	      
	      printf("trid[%d] inf-norm Aev-lev = %.16e\n",j,alpha);
	      printf("min eval = %.16e\n",beta);
#endif
	      bss_free(temp);
	      break;
	    }
	  beta = hmt_ddot(temp+1,r_k+1,k)/hmt_ddot(r_k+1,r_k+1,k);
	  /* beta = fddot(k,temp+1,1,r_k+1,1)/fddot(k,r_k+1,1,r_k+1,1); */
	  bss_free(temp);
	  rr=0.0;
	}
      else
	{rr+=1.0e-14;}
      if (j>10)
	{
#ifdef DEBUG      
	  printf("rr=%.16e\n",rr);
	  printf("beta=%.16e\n",beta);
#endif
	  break;
	}
    }
  fflush(stdout);
  
  if (k>1)
    {
      beta_old = dv[2];
#ifdef DEBUG      
      printf("begin seed vector calc\n");
#endif
      /* tri_calc_evec(dv_w,ev_w,k,beta_old,r_k_old); */
      r_k_old[0]=0.0;
      dvec_set(r_k_old+1,1.0,k);
      normalize(r_k_old+1,k);
      j=0;
      temp=trid_mult(dv_w,ev_w,k,r_k_old);
      alpha_old=0.0;
      for (i=1; i<=k; i++)
	{alpha_old = MAX(alpha,fabs(beta_old*r_k_old[i]-temp[i]));}
      bss_free(temp);
      /* printf("trid[%d] inf-norm Aev-lev = %.16e\n",j,alpha_old); */
      rr = 0.0;
      for (;;)
	{
	  j++;
	  if (trid_slv(dv_w,ev_w,k,beta_old-rr,r_k_old))
	    {
	      temp = trid_mult(dv_w,ev_w,k,r_k_old);
	      alpha_old=0.0;
	      for (i=1; i<=k; i++)
		{alpha_old = MAX(alpha,fabs(beta_old*r_k_old[i]-temp[i]));}

	      if ((alpha_old<1.0e-15)||(j))
		{
#ifdef DEBUG
		  printf("trid[%d] inf-norm Aev-lev = %.16e\n",j,alpha_old);
		  printf("min eval = %.16e\n",beta_old);
#endif
		  bss_free(temp);
		  break;
		}
	      beta_old = hmt_ddot(temp+1,r_k_old+1,k)/
		hmt_ddot(r_k_old+1,r_k_old+1,k);
	      /* beta_old = fddot(k,temp+1,1,r_k_old+1,1)/
		 fddot(k,r_k_old+1,1,r_k_old+1,1);  */
	      bss_free(temp);
	      rr = 0.0;
	    }
	  else
	    {rr+=1.0e-14;}
	  if (j>10)
	    {
#ifdef DEBUG      
	      printf("rr=%.16e\n",rr);
	      printf("beta_old=%.16e\n",beta_old);
#endif
	      break;
	    }
	}
    }
  fflush(stdout);

  dvec_zero(ev,n+1);
  dvec_zero(rhs,n+1);
  dvec_zero(new_rhs,2*n);


  for (j=1; j<=k; j++)
    {
      if (k==1)
	{Init_rhs_random(rhs,n);}
      else
	{hmt_daxpy(rhs,r_k_old[j],Q[j],n);}
      hmt_daxpy(new_rhs,r_k[j],Q[j],n);
      /* {fdaxpy(n,r_k_old[j],Q[j],1,rhs,1);} */
      /* fdaxpy(n,r_k[j]    ,Q[j],1,new_rhs ,1); */

      /*

      for (i=0; i<n; i++)
	{

	  rhs[i]  += r_k_old[j] * Q[j][i];
	  ev[i]   += r_k[j] * Q[j][i];
	  hmt_daxpy(rhs,r_k_old[j],Q[j],n);
	  hmt_daxpy(ev ,r_k[j],Q[j],n);



	  this was the good part of the for loop
	  rhs[i]  += r_k_old[j] * Q[j][i];
	  ev[i]   += r_k[j] * Q[j][i];

	}
      */
    }
  normalize(ev,n);
  normalize(new_rhs,n);
  normalize(rhs,n);

  /*check A evec_min = eval_min evec_min */
  matrix_mult(Ap_k,new_rhs,n,ia,ja,vals);
  alpha=0.0;
  for (i=0; i<n; i++)
    {alpha = MAX(alpha,fabs(beta*new_rhs[i]-Ap_k[i]));}
#ifdef DEBUG      
  printf("1 :: fiedler tqli inf-norm Aev-lev = %.16e\n",alpha);

  if (alpha>1.0e-01)
    {printf("poor fiedler vector - increase res_limit.\n");}
#endif

  /*check A evec_min = eval_min evec_min */
  matrix_mult(Ap_k,rhs,n,ia,ja,vals);
  alpha=0.0;
  for (i=0; i<n; i++)
    {alpha = MAX(alpha,fabs(beta_old*rhs[i]-Ap_k[i]));}

#ifdef DEBUG      
  printf("2 :: new rhs tqli inf-norm Aev-lev = %.16e\n",alpha);

  if (alpha>1.0e-01)
    {printf("poor seed for next round - increase res_lim.\n");}
#endif

  while (alpha<1.0e-07)
    {
      hmt_perturb(rhs,n);
      matrix_mult(Ap_k,rhs,n,ia,ja,vals);
      alpha=0.0;
      for (i=0; i<n; i++)
	{alpha = MAX(alpha,fabs(beta_old*rhs[i]-Ap_k[i]));}
    }
#ifdef DEBUG      
  printf("3 :: new rhs tqli inf-norm Aev-lev = %.16e\n",alpha);
#endif

  
  for (i=1;i<=k;i++)
    {perm_free(Q[i]);}
  perm_free(Q);

  /* return is in r_k_new so don't free */
  bss_free(r_k);
  bss_free(r_k_old);
  bss_free(r_k_new);
  bss_free(p_k);
  bss_free(x_k);
  bss_free(Ap_k);

#ifdef DEBUG
  printf("perm diff = %d\n",(perm_calls()-perm_frees())-p);
  printf("bss  diff = %d\n",(bss_calls()-bss_frees())-b);
#endif

  return(new_rhs);
}



/********************************sparse_matrix.c*******************************
Function: 

Input : 
Output: 
Return: 
Description:  
*********************************sparse_matrix.c******************************/
struct tree_node *
new_tree_node()
{
  struct tree_node *new;
  int size;

  
  size = sizeof(struct tree_node);

#ifdef DEBUG
  if (!size)
    {error_msg_fatal("new_tree_node() :: structure is zero bytes in size!\n");}
#endif


  new = (struct tree_node *) bss_malloc(size);

  /* initialize to all zeros */
  /*
  if (!(size%sizeof(double)))
    {dvec_zero((double *)new,size/sizeof(double));}
  else if (!(size%INT_LEN))
    {ivec_zero((INT *)new,size/INT_LEN);}
  */
  if (!(size%sizeof(char)))
    {bzero((char *)new,size/sizeof(char));}
  else
    {error_msg_fatal("new_tree_node() :: can't initialize tree_node!\n");}

  return(new);
}


/********************************sparse_matrix.c*******************************
Function: 

Input : 
Output: 
Return: 
Description:  
*********************************sparse_matrix.c******************************/
int
log_2_floor(int n)
{
  int i=1, depth=0;

  if (n<1) {error_msg_fatal("log_2_floor() :: n>=1!!!");}
  while (i<=n)
    {depth++; i<<=1;}

  return(--depth);
}





/********************************sparse_matrix.c*******************************
Function: 

Input : 
Output: 
Return: 
Description:  
*********************************sparse_matrix.c******************************/
/*
matrix_id_ptr
extract_sub_matrix(int *m,int n)
{
;
}
*/


/********************************sparse_matrix.c*******************************
Function: 

Input : 
Output: 
Return: 
Description:  
*********************************sparse_matrix.c******************************/
void
compress_adj_list(void)
{
  int  *lens;
  int **rows;
  double **vals;
  int i, j, n, rl;
  int *r;
  double *v;
  int elmt;
  adj_node_ptr base, offset;


  n = active_matrix->lda;
  lens = (int *)     bss_malloc((n+1)*sizeof(int));
  rows = (int **)    bss_malloc((n+1)*sizeof(int *));
  vals = (double **) bss_malloc((n+1)*sizeof(double *));

  base = active_matrix->adj_list;
  for (i=0;i<n;i++)
    {
      lens[i] = rl = base->elmt;
      if (rl<2)
	{printf("compress_adj_list() :: one or fewer entries?");}
      rows[i] = r  = (int *)    bss_malloc((rl+1)*sizeof(int));
      vals[i] = v  = (double *) bss_malloc(rl*sizeof(double));
      offset = (base++)->next;
      for (j=0;j<rl;j++)
	{
	  *r++ = elmt = offset->elmt;
	  *v++ = offset->val;
	  offset = offset->next;
	  /*
	  if (elmt==i)
	    {*v++ = 1.0*(rl-1);}
	  else
  	    {*v++ = -1.0;}
	  */
	}
      *r = -1;
    }
  lens[n] = -1;
  rows[n] = NULL;
  vals[n] = NULL;
  active_matrix->row_lens    = lens;
  active_matrix->row_indices = rows;
  active_matrix->row_vals    = vals;
  

#ifdef DEBUG_0
  printf("\n");
  for (i=0;i<n;i++)
    {
      r = rows[i];
      v = vals[i];
      for (j=0;j<lens[i];j++)
	{printf("%3d ",*r++ + 1);}
      printf("\n");
      for (j=0;j<lens[i];j++)
	{printf("%3.0f ",*v++);}
      printf("\n");
    }
#endif
}

/********************************sparse_matrix.c*******************************
Function: SMI_init

Input : 
Output: 
Return: 
Description:  
*********************************sparse_matrix.c******************************/
void
SMI_init(void)
{
int i;
void SMI_free_matrix(int);


/* if this is not the first call */
if (num_mat)
  {
     for (i=1;i<=MAX_MAT;i++)
       {SMI_free_matrix(i);}
   }

/* initalize */
if (num_mat)
  {error_msg_fatal("SMI_init has num_mat error");}
else
  {
    for (i=0;i<MAX_MAT;i++)
      {matrix_list[i] = NULL;}
  }
}



/********************************sparse_matrix.c*******************************
Function: SMI_free_matrix

Input : 
Output: 
Return: 
Description:  
*********************************sparse_matrix.c******************************/
void
SMI_free_matrix(int mat)
{
  matrix_id_ptr *id;
  int i, lda;

  if ((mat>=1)&&(mat<=MAX_MAT))
    {
      /* adjust 1...k to 0...k-1 */
      mat--;

      /* need to free adj_list, permutation, row, and columns!!! later!!! */
      /* possibly bad matrix id */
      id = &matrix_list[mat];
      if (*id)
	{
	  num_mat--;
	  lda = (*id)->lda;
	  free_adj_list((*id)->adj_list,lda);

	  for (i=0;i<lda;i++)
	    {
	      bss_free((*id)->row_indices[i]);
	      bss_free((*id)->row_vals[i]);
	    }

	  bss_free((*id)->row_lens);
	  bss_free((*id)->row_indices);
	  bss_free((*id)->row_vals);

	  bss_free(*id);
	  *id = NULL;
	}
      else
	{return;}
    }
  else
    {error_msg_fatal("trying to return an invalid matrix id?!?");}
}



/********************************sparse_matrix.c*******************************
Function: SMI_get_matrix

Input : 
Output: 
Return: 
Description:  
*********************************sparse_matrix.c******************************/
matrix_id_ptr
SMI_get_matrix(int *rt)
{
*rt=0;

/* if there are vacancies */
if (num_mat<MAX_MAT)
  {
     /* find the first - linear search */
     while ((*rt<MAX_MAT)&&(matrix_list[*rt] != NULL))
       {(*rt)++;}

     if (*rt<MAX_MAT)
       {
         ++num_mat;
         return(matrix_list[(*rt)++]=(matrix_id_ptr)bss_malloc(sizeof(matrix_id)));
	}
   }
*rt=ERROR; 
return(NULL);
}


/* prototypes */
static int load_file(FILE *);
static int pass_one(int *,FILE *);
static int f_pos(FILE *, int);
static int f_pos(FILE *, int);

/********************************sparse_matrix.c*******************************
Function: SMI_read_matlab()

Input : 
Output: 
Return: 
Description:  
*********************************sparse_matrix.c******************************/
int
SMI_read_matlab(char *file_name)
{
matrix_id_ptr id, SMI_get_matrix(int *);
int lda = 0;
int rt = 0;
int nnz = 0;
int i;
int *iptr;
double *dptr;
FILE *ifp;


  /* new template and add to list */
  if ((id = SMI_get_matrix(&rt)) == NULL)
    {return(ERROR);}

  strcpy (id->ty,"Sparse X");
  strcpy (id->desc,file_name);
    
  /* need file pointer to sparse matrix matlab formatted file */
  if ((ifp = fopen(file_name,"r")) == NULL)
    {error_msg_fatal("Can't open input file!!!");}
  
  lda = load_file(ifp);
  id->lda = lda;
  id->sda = sda;
 
  id->adj_list = init_adj_list(lda);

  lda =  *(pairs+1);
  iptr = (pairs+2);
  dptr = vals;
  for (i=0; i<lda; i++)
    {
      /* strip out zero - true sparse matrix only */
      if (*dptr != 0.0)
	{
	  insert_adj_list(id->adj_list,*iptr-1,*(iptr+1)-1,*dptr);
	  iptr+=2;
	  dptr++;
	  nnz++;
	}
    }

  bss_free(pairs);
  bss_free(vals);

  if (fclose(ifp))
    {error_msg_fatal("Can't close input file!!!");}

  id->permutation = NULL;
  id->row_vals = NULL;
  id->col_vals = NULL;
  id->nnz = nnz;
  return(rt);
}



/********************************sparse_matrix.c*******************************
Function: SMI_fprint_id_hb_no_vals(matrix,FULL,"k6112.hb")

Input : 
Output: 
Return: 
Description:  
*********************************sparse_matrix.c******************************/
int
SMI_fprint_id_hb_no_vals(int mat, int type,char *name)
{
matrix_id_ptr id;
char buf[81];
FILE *tfp, *ofp;
int len;

if ((mat<0)||(mat>MAX_MAT)||((id = matrix_list[mat-1])==NULL))
  {error_msg_fatal("Invalid matrix handle!!!");}


/* write strings to tmp file */
/*
if ((tfp = fopen("temp_file","w+")) == NULL)
  {error_msg_fatal("Can't open temp_file!!!");} 
*/

if ((tfp = tmpfile()) == NULL)
  {error_msg_fatal("Can't open temp_file!!!");} 

fprintf(tfp," (%d) ",mat);  
fprintf(tfp," %dx%d ",id->lda,id->lda);

rewind(tfp);
if (fgets(buf,80,tfp)==NULL)
  {error_msg_fatal("Can't read from temp_file!!!");} 
if (fclose(tfp))
  {error_msg_fatal("Can't close temp_file!!!");} 

if ((int)strlen(buf)>80)
  {printf("Harwell-Boeing Title > 80 Char Long!!!");}

if (!(ofp=fopen(name,"w")))
  {error_msg_fatal("Can't open hb file!!!");} 

fprintf(ofp,"%80s\n",buf);

/*
if (id->permutation == NULL)
  {printf("\n%s","Per: NULL");}
if (id->adj_list == NULL)
  {printf("\n%s","Adj: NULL");}
if (id->row_vals == NULL)
  {printf("\n%s","RV : NULL");}
if (id->col_vals == NULL)
  {printf("\n%s","CV : NULL");}
*/


/* remaining three lines of header determinned here!!! */
fprint_adj_list_hb_no_vals(ofp,id->adj_list, id->lda, id->sda, type);

if (fclose(ofp))
  {error_msg_fatal("Can't close hb file!!!");} 

return(mat);
}





/********************************sparse_matrix.c*******************************
Function: SMI_print_id_hb

Input : 
Output: 
Return: 
Description:  
*********************************sparse_matrix.c******************************/
int
SMI_print_id_hb(int mat, int type)
{
matrix_id_ptr id;
char buf[81];
FILE *tfp;
int len;

if ((mat<0)||(mat>MAX_MAT)||((id = matrix_list[mat-1])==NULL))
  {error_msg_fatal("Invalid matrix handle!!!");}


/* write strings to tmp file */
/*
if ((tfp = fopen("temp_file","w+")) == NULL)
  {error_msg_fatal("Can't open temp_file!!!");} 
*/

if ((tfp = tmpfile()) == NULL)
  {error_msg_fatal("Can't open temp_file!!!");} 

fprintf(tfp," (%d) ",mat);  
fprintf(tfp," %dx%d ",id->lda,id->lda);
fprintf(tfp,"%s<-",id->ty);
fprintf(tfp,"%s",id->desc);

rewind(tfp);
if (fgets(buf,80,tfp)==NULL)
  {error_msg_fatal("Can't read from temp_file!!!");} 
if (fclose(tfp))
  {error_msg_fatal("Can't close temp_file!!!");} 

if ((int)strlen(buf)>80)
  {printf("Harwell-Boeing Title > 80 Char Long!!!");}
printf("%80s\n",buf);

/*
if (id->permutation == NULL)
  {printf("\n%s","Per: NULL");}
if (id->adj_list == NULL)
  {printf("\n%s","Adj: NULL");}
if (id->row_vals == NULL)
  {printf("\n%s","RV : NULL");}
if (id->col_vals == NULL)
  {printf("\n%s","CV : NULL");}
*/


/* remaining three lines of header determinned here!!! */
print_adj_list_hb(id->adj_list, id->lda, id->sda, type);

return(mat);
}



/********************************driver.c**************************************
Function load_file() 

Input : pointers to i/o space.
Output: allocates space and loads data.
Return: none.
Description: allocates space for the program and reads in data from
input file specified. error checking done.
*********************************driver.c*************************************/
static
int 
load_file(FILE *ifp)
{
int i, j, n = -1, count = 0, vc = 0;
char num[MAX_LEN], *ptr;
double v;


/* first pass through input                     */
/* count the number of (i,j) pair minus (i,i)'s */
/* n = get LDA of matrix                        */
n=pass_one(&count,ifp);


/* allocate space                                     */
/* n = LDA -> space to allocate for permutation array */
/* total space for input is 2*count +2                */
if ((n <= 0) || (n > MAX_LDA))
  {error_msg_fatal("LDA <= 0 OR LARGER THAN MAX_LDA !!!");}
if ((count <= 0) || (count > n*n))
  {error_msg_fatal("(i,j) entries <= 0 OR LARGER THAN LDA^2 !!!");}

if ((pairs = (int *) bss_malloc((2*count+2)*INT_LEN)) == NULL)
  {error_msg_fatal("Could not allocate space for pairs!?!");}
if ((vals = (double *) bss_malloc(count*sizeof(double))) == NULL)
  {error_msg_fatal("Could not allocate space for values!?!");}


/* load n and 2*count into first two positions */
*(pairs+0)=n;
*(pairs+1)=count;

/* skip in ... */
count = 2;

/* start over again ... */
rewind(ifp);

/* second pass         */
/* read in count pairs */
/* transposes X so we read in row- wise */
while (fgets(num,MAX_LEN,ifp) != NULL)
  {
     /* get first of (i,j) pair */
     if ((ptr=(char *) strtok(num,DELIM)) != NULL)
       {
	 i = atoi(ptr);
	 /*	 j = atoi(ptr); */
    
	 /* second of pair */
	 if ((ptr=(char *) strtok(NULL,DELIM)) != NULL)
	   {
	     /*	     i = atoi(ptr); */
	     j = atoi(ptr);
	     
	     /* value */
	     if ((ptr=(char *) strtok(NULL,DELIM)) != NULL)
	       {v = atof(ptr);}
	   }
       }
   
     /* add to pair and val list */
     *(pairs+count)=i;
     count++;
     *(pairs+count)=j;
     count++;
     *(vals+vc) = v;
     vc++;
   }


if (vc!=*(pairs+1))
  {printf("vc off!!!");}

return(n);
}




/********************************driver.c**************************************
Function pass_one()

Input : address of count.
Output: number of (i,j) pairs
Return: LDA of matrix
Description: processes input to find LDA and the number of (i,j |i!=j)
pairs.
*********************************driver.c*************************************/
static
int 
pass_one(int *ct, FILE *ifp)
{

int i, j, n=0, m=0, count=0;
char num[MAX_LEN], *ptr;



/* first pass through input                     */
/* count the number of (i,j) pair minus (i,i)'s */
while (fgets(num,MAX_LEN,ifp) != NULL)
  {
     /* get first of (i,j) pair */
     if ((ptr=(char *) strtok(num,DELIM)) != NULL)
       {
          i = atoi(ptr);
          if ((i<0)||(i>MAX_LDA))
            {error_msg_fatal("BAD i !!!");}
          if (i>n)
	    {n=i;}
    
          /* second of pair */
          if ((ptr=(char *) strtok(NULL,DELIM)) != NULL)
	    {
               j = atoi(ptr);
               if ((j<0)||(j>MAX_LDA))
                 {error_msg_fatal("BAD j !!!");}
               if (j>m)
         	 {m=j;}
	     }
          else
	    {error_msg_fatal("j INPUT EXPECTED !!!");}
	}
     else
       {error_msg_fatal("i INPUT EXPECTED !!!");}

    
     count++;
   }

/*
if (n!=m)
  {printf("matlab file is not square!");}
*/

/* LATER: a real hack but I'll fix it later */
sda = m;

/* return #(i,j) pairs and LDA */
*ct=count;
return(n);
}

/********************************sparse_matrix.c*******************************
Function: SMI_read_hb

Input : 
Output: 
Return: 
Description:  
*********************************sparse_matrix.c******************************/
int
SMI_read_hb(char *file_name, int hack)
{
matrix_id_ptr id, SMI_get_matrix(int *);
char input[80], *c_ptr;
int tot_lines = 0, num_ptr_lines, num_vert_lines, num_val_lines, num_vec_lines;
int lda = 0, sda, nnz=0, nvnz=0;
int rt = 0;
int i,j,k;
int *ptrs, *cols;
double *vals;
FILE *ifp;
int *tmp;


/* new template and add to list */
if ((id = SMI_get_matrix(&rt)) == NULL)
  {error_msg_fatal("No matrix handles left!!! Increase MAX_MAT.");}


/* need file pointer */
if ((ifp = fopen(file_name,"r")) == NULL)
  {error_msg_fatal("Can't open input file!!!");}

/* reading the header should be another function but ... started bnf */
/* the header of a HB ascii file is four lines	*/
/* the first contains the name of the HB file	*/
if (fgets(input,MAX_LEN,ifp) != NULL)
  {
    if ((c_ptr=(char *) strtok(input,DELIM)) != NULL)
      {strncpy(id->desc,c_ptr,80);}
    else
      {error_msg_fatal("First Line of HB file is blank?!?");}
     

     /* remove cr/lf */
     /*
     sda = strcspn(input,"\0");
     memcpy(id->desc,input,--sda);
     */
   }
else
  {rt=ERROR; printf("\nBAD FIRST LINE HB FILE - FILE EMPTY?");}

/* the second contains number of lines of following sections */
if (fgets(input,MAX_LEN,ifp) != NULL)
  {
     if ((c_ptr=(char *) strtok(input,DELIM)) != NULL)
       {tot_lines = atoi(c_ptr);}
     if ((c_ptr=(char *) strtok(NULL,DELIM)) != NULL)
       {num_ptr_lines = atoi(c_ptr);}
     if ((c_ptr=(char *) strtok(NULL,DELIM)) != NULL)
       {num_vert_lines = atoi(c_ptr);}
     if ((c_ptr=(char *) strtok(NULL,DELIM)) != NULL)
       {num_val_lines = atoi(c_ptr);}
     if ((c_ptr=(char *) strtok(NULL,DELIM)) != NULL)
       {num_vec_lines = atoi(c_ptr);}
     if ((!tot_lines)||
         (tot_lines!=num_ptr_lines+num_vert_lines+num_val_lines+num_vec_lines))
       {rt=ERROR; printf("\nBAD SECOND LINE HB FILE - LINES DON'T ADD UP!");}
   }
else
  {rt=ERROR; printf("\nBAD SECOND LINE HB FILE - ONLY ONE LINE IN HB FILE?");}

/* third line contains matrix dimensions and number of	*/
/* non-zeros in the matrix upper half -> symmetric!!!	*/
/* must match types! what about non_square matrix?	*/
if (fgets(input,MAX_LEN,ifp) != NULL)
  {
     if ((c_ptr=(char *) strtok(input,DELIM)) != NULL)
       {strcpy(id->ty,c_ptr);}
     if ((c_ptr=(char *) strtok(NULL,DELIM)) != NULL)
       {lda = atoi(c_ptr);}
     if ((c_ptr=(char *) strtok(NULL,DELIM)) != NULL)
       {sda = atoi(c_ptr);}
     if ((c_ptr=(char *) strtok(NULL,DELIM)) != NULL)
       {nnz = atoi(c_ptr);}
     if ((c_ptr=(char *) strtok(NULL,DELIM)) != NULL)
       {nvnz = atoi(c_ptr);}
     /*
     if (lda!=sda)
       {printf("HB file does does not contain a square matrix!");}
     */
     if ((!lda)||(!sda)||(!nnz))
       {rt=ERROR; printf("\nBAD THIRD LINE HB FILE");}
   }
else
  {rt=ERROR; printf("\nBAD THIRD LINE HB FILE - ONLY TWO LINES IN HB FILE?");}

/* fourth line contains the format statements but who cares?	*/
/* pointer and vert values must be integers -> ignore format	*/
/* matrix and vector values can be either int or double	however	*/
/* we're going to make them both double and let the user know.	*/
if (fgets(input,MAX_LEN,ifp) != NULL)
  {
     /* so ignore the first two */
     if ((c_ptr=(char *) strtok(input,DELIM)) != NULL)
     if ((c_ptr=(char *) strtok(NULL,DELIM)) != NULL)

     /* let user know about matrix coercion */
     if (num_val_lines > 0)
       {
          if ((c_ptr=(char *) strtok(NULL,DELIM)) != NULL)
	    {
               if ((strchr(c_ptr,'i')!=NULL)||(strchr(c_ptr,'I')!=NULL))
		 {printf("\nCoercing to doubles: %s",c_ptr);}
	     }
	}

     /* let user know about matrix coercion */
     if (num_vec_lines > 0)
       {
          if ((c_ptr=(char *) strtok(NULL,DELIM)) != NULL)
	    {
               if ((strchr(c_ptr, 'i')!=NULL)||(strchr(c_ptr,'I')!=NULL))
		 {printf("\nCoercing to doubles: %s",c_ptr);}
            }
	}
   }
else
  {
     rt = ERROR; 
     printf("\nBAD FOURTH LINE HB FILE - ONLY THREE LINES IN HB FILE?");
   }


/* time to make the donuts */
/* must be hungry I mean time to read the matrix in */
if (f_pos(ifp,HEADER_LEN))
  {
    /* space to read lists in to temp storage */
    ptrs = (int *) bss_malloc((lda+1)*sizeof(int));
    cols = (int *) bss_malloc(nnz*sizeof(int));
    vals = (double *) bss_malloc(nnz*sizeof(double));
    
    if ((ptrs==NULL)||(cols==NULL)||(vals==NULL))
      {error_msg_fatal("not enough space in read_hb!!!");}

    /* pointers */
    j=0;
    for (i=0; i<num_ptr_lines; i++)
      {
	fgets(input,MAX_LEN,ifp);

        *(ptrs+j) = atoi(strtok(input,DELIM));
        j++;
  
	while ((c_ptr=(char *) strtok(NULL,DELIM))!=NULL)
         {*(ptrs+j) = atoi(c_ptr); j++;}
      }

    if (j!=(lda+1))
      {error_msg_fatal("Actual Pointer Count differs from lda+1!!!");}

    /* vertices */
    j=0;
    for (i=0; i<num_vert_lines; i++)
      {
	fgets(input,MAX_LEN,ifp);

        *(cols+j) = atoi(strtok(input,DELIM));
        j++;
  
	while ((c_ptr=(char *) strtok(NULL,DELIM))!=NULL)
         {*(cols+j) = atoi(c_ptr); j++;}
      }


    if (j!=nnz)
      {error_msg_fatal("Vertice Count differs from nnz!!!");}

    /* values */
    j=0;

    for (i=0; i<num_val_lines; i++)
      {
	fgets(input,MAX_LEN,ifp);

        *(vals+j) = atof(strtok(input,DELIM));
        j++;
  
	while ((c_ptr=(char *) strtok(NULL,DELIM))!=NULL)
         {*(vals+j) = atof(c_ptr); j++;}
      }

    if (j&&(j!=nnz))
      {error_msg_fatal("Value Count differs from nnz!!!");}

    if (!j)
      {
	printf("No Values");
	for (i=0;i<nnz;i++)
	  {*(vals+j) = 0.0;}
      }
    
    /* initialize adj list */
    id->adj_list = init_adj_list(lda);

    /* place i,j pairs into adj list */
    /* symmetric case but don't add diags twice */
    if (((strchr(id->ty,'S')!=NULL)||(strchr(id->ty,'s')!=NULL))&&(hack==SYMM))
      {
	k=0;
	for (i=0; i<lda; i++)
	  {
	    for (j=0; j<(*(ptrs+i+1)-*(ptrs+i)); j++)
	      {
		insert_adj_list(id->adj_list,i,*(cols+k)-1,*(vals+k)); 
		if (i!=*(cols+k)-1)
		  {insert_adj_list(id->adj_list,*(cols+k)-1,i,*(vals+k));}
		k++;
	      }
	  }
      }
    else
      {
	k=0;
	for (i=0; i<lda; i++)
	  {
	    for (j=0; j<(*(ptrs+i+1)-*(ptrs+i)); j++)
	      {
		insert_adj_list(id->adj_list,i,*(cols+k)-1,*(vals+k)); 
		k++;
	      }
	  }
      }


    if (k!=nnz)
      {error_msg_fatal("didn't place nnz entries in adj_list");}

    /* free temp storage */
    bss_free(ptrs);
    bss_free(cols);
    bss_free(vals);
   }

/* the file is too short or number of line entries bad */
else
  {rt=ERROR; printf("\nCan't Position Someone");}

if (fclose(ifp))
  {error_msg_fatal("Can't close input file!!!");}

id->lda = lda;
id->sda = sda;
id->nnz = nnz;
id->permutation = NULL;
id->row_vals = NULL;
id->col_vals = NULL;

return(rt);
}



/********************************sparse_matrix.c*******************************
Function: f_pos()

Input : 
Output: 
Return: 
Description:  
*********************************sparse_matrix.c******************************/
static 
int
f_pos(FILE *fp, int num_lines)
{
int i;
char input[MAX_LEN];

rewind(fp);

for (i=0; i<num_lines; i++)
  {
     if (fgets(input,MAX_LEN,fp) == NULL)
       {return(FAIL);}
   }

return(PASS);
}



/********************************sparse_matrix.c*******************************
Function: SMI_print_id

Input : 
Output: 
Return: 
Description:  
*********************************sparse_matrix.c******************************/
int
SMI_print_id_matlab(int mat)
{
matrix_id_ptr id;

if ((mat<0)||(mat>MAX_MAT)||((id = matrix_list[mat-1])==NULL))
  {error_msg_fatal("Invalid matrix handle!!!");}

#ifdef INFO
printf("\n\nMATRIX ID: %d",mat);
printf("\nLDA: %d",id->lda);
printf("\nType: %s",id->ty);
printf("\nDesc:\n%s",id->desc);
if (id->permutation == NULL)
  {printf("\n%s","Per: NULL");}
if (id->adj_list == NULL)
  {printf("\n%s","Adj: NULL");}
if (id->row_vals == NULL)
  {printf("\n%s","RV : NULL");}
if (id->col_vals == NULL)
  {printf("\n%s","CV : NULL");}
#endif


print_adj_list_matlab(id->adj_list,id->lda);

return(mat);
}



/********************************sparse_matrix.c*******************************
Function: SMI_make_symm

Input : 
Output: 
Return: 
Description:  
*********************************sparse_matrix.c******************************/
int
SMI_make_symm(int mat)
{
matrix_id_ptr id;
int len;

if ((mat<0)||(mat>MAX_MAT)||((id = matrix_list[mat-1])==NULL))
  {return(ERROR);}


/* remaining three lines of header determinned here!!! */
if (!is_adj_list_symm(id->adj_list, id->lda))
  { 
    if (make_adj_list_symm(id->adj_list, id->lda)==PASS)
      {return(SYMM);}
    else
      {printf("found to be symmetric but couldn't symmatrize?!?");}
  }

return(ERROR);
}




/********************************sparse_matrix.c*******************************
Function: SMI_get_dim

Input : 
Output: 
Return: 
Description:  
*********************************sparse_matrix.c******************************/
int
SMI_get_dim(int mat, int *lda, int *sda)
{
  matrix_id_ptr id;


  if ((mat<0)||(mat>MAX_MAT)||((id = matrix_list[mat-1])==NULL))
    {return(ERROR);}


  *lda = id->lda;
  *sda = id->sda;

  return(PASS);
}



/********************************sparse_matrix.c*******************************
Function: SMI_extract_row

Input : 
Output: 
Return: 
Description:  
*********************************sparse_matrix.c******************************/
adj_node_ptr
SMI_extract_row(int mat, int row)
{
  matrix_id_ptr id;


  if ((mat<0)||(mat>MAX_MAT)||((id = matrix_list[mat-1])==NULL))
    {return(NULL);}
  if ((row<0)||(row>=id->lda))
    {return(NULL);}

  return(id->adj_list+row);
}



/********************************sparse_matrix.c*******************************
Function: SMI_extract_matrix

Input : 
Output: 
Return: 
Description:  
*********************************sparse_matrix.c******************************/
matrix_id_ptr
SMI_extract_matrix(int mat)
{
  matrix_id_ptr id;


  if ((mat<0)||(mat>MAX_MAT)||((id = matrix_list[mat-1])==NULL))
    {return(NULL);}
 
  return(id);
}




/********************************sparse_matrix.c*******************************
Function: SMI_get_nnz

Input : 
Output: 
Return: 
Description:  
*********************************sparse_matrix.c******************************/
int
SMI_get_nnz(int mat, int *nnz)
{
  matrix_id_ptr id;
  int i, sum=0;

  if ((mat<0)||(mat>MAX_MAT)||((id = matrix_list[mat-1])==NULL))
    {return(ERROR);}

  if (id->nnz>0)
    {
      *nnz=id->nnz;
      return(PASS);
    }

  /* if we get here id->adj_list == NULL */
  error_msg_fatal("nnz negative?!? did you initialize?!?");
  return(FAIL);
}



/********************************sparse_matrix.c*******************************
Function: SMI_get_adj

Input : 
Output: 
Return: 
Description:  
*********************************sparse_matrix.c******************************/
adj_node_ptr
SMI_get_adj(int mat)
{
  matrix_id_ptr id;
  int i, sum=0;

  if ((mat<0)||(mat>MAX_MAT)||((id = matrix_list[mat-1])==NULL))
    {return(NULL);}

  return(id->adj_list);
}




/********************************sparse_matrix.c*******************************
Function: SMI_read_crs()

Input : 
Output: 
Return: 
Description:  
*********************************sparse_matrix.c******************************/
int
SMI_read_crs(char *file_name, int gs, int hack)
{
matrix_id_ptr id, SMI_get_matrix(int *);
char input[80], *c_ptr;
int tot_lines = 0, num_ptr_lines, num_vert_lines, num_val_lines, num_vec_lines;
int lda = 0, sda, nnz=0, nvnz=0;
int rt = 0;
int i,j,k;
int *ptrs, *cols, *iptr, **i_next, *c_next;
double *vals;
FILE *ifp;
int *tmp;
int row_len;
int *row;
double *row_vals, *dptr, **v_next;

/* new template and add to list */
if ((id = SMI_get_matrix(&rt)) == NULL)
  {return(ERROR);}

/* need file pointer */
if ((ifp = fopen(file_name,"r")) == NULL)
  {error_msg_fatal("read_crs_hb :: can't open input file!!!");}

/* reading the header should be another function but ... started bnf */
/* the header of a HB ascii file is four lines	*/
/* the first contains the name of the HB file	*/
if (fgets(input,MAX_LEN,ifp) != NULL)
  {
    if ((c_ptr=(char *) strtok(input,DELIM)) != NULL)
      {strncpy(id->desc,c_ptr,80);}
    else
      {error_msg_fatal("First Line of HB file is blank?!?");}
     

     /* remove cr/lf */
     /*
     sda = strcspn(input,"\0");
     memcpy(id->desc,input,--sda);
     */
   }
else
  {rt=ERROR; printf("\nBAD FIRST LINE HB FILE - FILE EMPTY?");}

/* the second contains number of lines of following sections */
if (fgets(input,MAX_LEN,ifp) != NULL)
  {
     if ((c_ptr=(char *) strtok(input,DELIM)) != NULL)
       {tot_lines = atoi(c_ptr);}
     if ((c_ptr=(char *) strtok(NULL,DELIM)) != NULL)
       {num_ptr_lines = atoi(c_ptr);}
     if ((c_ptr=(char *) strtok(NULL,DELIM)) != NULL)
       {num_vert_lines = atoi(c_ptr);}
     if ((c_ptr=(char *) strtok(NULL,DELIM)) != NULL)
       {num_val_lines = atoi(c_ptr);}
     if ((c_ptr=(char *) strtok(NULL,DELIM)) != NULL)
       {num_vec_lines = atoi(c_ptr);}
     if ((!tot_lines)||
         (tot_lines!=num_ptr_lines+num_vert_lines+num_val_lines+num_vec_lines))
       {rt=ERROR; printf("\nBAD SECOND LINE HB FILE - LINES DON'T ADD UP!");}
   }
else
  {rt=ERROR; printf("\nBAD SECOND LINE HB FILE - ONLY ONE LINE IN HB FILE?");}

/* third line contains matrix dimensions and number of	*/
/* non-zeros in the matrix upper half -> symmetric!!!	*/
/* must match types! what about non_square matrix?	*/
if (fgets(input,MAX_LEN,ifp) != NULL)
  {
     if ((c_ptr=(char *) strtok(input,DELIM)) != NULL)
       {strcpy(id->ty,c_ptr);}
     if ((c_ptr=(char *) strtok(NULL,DELIM)) != NULL)
       {lda = atoi(c_ptr);}
     if ((c_ptr=(char *) strtok(NULL,DELIM)) != NULL)
       {sda = atoi(c_ptr);}
     if ((c_ptr=(char *) strtok(NULL,DELIM)) != NULL)
       {nnz = atoi(c_ptr);}
     if ((c_ptr=(char *) strtok(NULL,DELIM)) != NULL)
       {nvnz = atoi(c_ptr);}
     /*
     if (lda!=sda)
       {printf("HB file does does not contain a square matrix!");}
     */
     if ((!lda)||(!sda)||(!nnz))
       {rt=ERROR; printf("\nBAD THIRD LINE HB FILE");}
   }
else
  {rt=ERROR; printf("\nBAD THIRD LINE HB FILE - ONLY TWO LINES IN HB FILE?");}

/* fourth line contains the format statements but who cares?	*/
/* pointer and vert values must be integers -> ignore format	*/
/* matrix and vector values can be either int or double	however	*/
/* we're going to make them both double and let the user know.	*/
if (fgets(input,MAX_LEN,ifp) != NULL)
  {
     /* so ignore the first two */
     if ((c_ptr=(char *) strtok(input,DELIM)) != NULL)
     if ((c_ptr=(char *) strtok(NULL,DELIM)) != NULL)

     /* let user know about matrix coercion */
     if (num_val_lines > 0)
       {
          if ((c_ptr=(char *) strtok(NULL,DELIM)) != NULL)
	    {
               if ((strchr(c_ptr,'i')!=NULL)||(strchr(c_ptr,'I')!=NULL))
		 {printf("\nCoercing to doubles: %s",c_ptr);}
	     }
	}

     /* let user know about matrix coercion */
     if (num_vec_lines > 0)
       {
          if ((c_ptr=(char *) strtok(NULL,DELIM)) != NULL)
	    {
               if ((strchr(c_ptr, 'i')!=NULL)||(strchr(c_ptr,'I')!=NULL))
		 {printf("\nCoercing to doubles: %s",c_ptr);}
            }
	}
   }
else
  {
     rt = ERROR; 
     printf("\nBAD FOURTH LINE HB FILE - ONLY THREE LINES IN HB FILE?");
   }


/* time to make the donuts */
/* must be hungry I mean time to read the matrix in */
if (f_pos(ifp,HEADER_LEN))
  {
    /* space to read lists in to temp storage */
    ptrs = (int *) bss_malloc((lda+1)*sizeof(int));
    cols = (int *) bss_malloc(nnz*sizeof(int));
    vals = (double *) bss_malloc(nnz*sizeof(double));
    
    if ((ptrs==NULL)||(cols==NULL)||(vals==NULL))
      {error_msg_fatal("not enough space in read_hb!!!");}

    /* pointers */
    j=0;
    for (i=0; i<num_ptr_lines; i++)
      {
	fgets(input,MAX_LEN,ifp);

        *(ptrs+j) = atoi(strtok(input,DELIM));
        j++;
  
	while ((c_ptr=(char *) strtok(NULL,DELIM))!=NULL)
         {*(ptrs+j) = atoi(c_ptr); j++;}
      }

    if (j!=(lda+1))
      {error_msg_fatal("Actual Pointer Count differs from lda+1!!!");}

    /* vertices */
    j=0;
    for (i=0; i<num_vert_lines; i++)
      {
	fgets(input,MAX_LEN,ifp);

        *(cols+j) = atoi(strtok(input,DELIM));
        j++;
  
	while ((c_ptr=(char *) strtok(NULL,DELIM))!=NULL)
         {*(cols+j) = atoi(c_ptr); j++;}
      }


    if (j!=nnz)
      {error_msg_fatal("Vertice Count differs from nnz!!!");}

    /* values */
    j=0;
    for (i=0; i<num_val_lines; i++)
      {
	fgets(input,MAX_LEN,ifp);

        *(vals+j) = atof(strtok(input,DELIM));
        j++;
  
	while ((c_ptr=(char *) strtok(NULL,DELIM))!=NULL)
         {*(vals+j) = atof(c_ptr); j++;}
      }

    if (j!=nnz)
      {error_msg_fatal("Value Count differs from nnz!!!");}
    
    /* place in crs format */
    /* symmetric case we have to add all but diags twice ...   */
    /* probably have to go through adj_list to accomplish this */
    if (((strchr(id->ty,'S')!=NULL)||(strchr(id->ty,'s')!=NULL))&&(hack==SYMM))
      {error_msg_fatal("ummm ... we have yet to write this routine!");}
    else
      {
	
	k=0;
	i_next = id->col_indices = (int **) bss_malloc(sizeof(int *)*lda);
	v_next = id->col_vals = (double **) bss_malloc(sizeof(double *)*lda);
	c_next = id->col_lens = (int *) bss_malloc(INT_LEN*lda);
	for (i=0; i<lda; i++)
	  {
	    c_next[i] = row_len = (*(ptrs+i+1)-*(ptrs+i));
	    if (row_len < (SPARSE*sda))
	      {
		i_next[i] = iptr = (int *) bss_malloc(INT_LEN*row_len); 
		v_next[i] = dptr = (double *) bss_malloc(sizeof(double)*row_len); 
		for (j=0; j<row_len; j++)
		  {
		    *iptr = *(cols+k)-1;
		    iptr++;
		    *dptr = *(vals+k);
		    dptr++;
		    k++;
		  }
	      }
	    else
	      {
		i_next[i] = iptr = NULL;
		v_next[i] = dptr = (double *) bss_malloc(sizeof(double)*sda); 
		for (j=0; j<sda; j++)
		  {
		    if ((cols[k]-1) == j)
		      {
			*dptr = *(vals+k);
			k++;
		      }
		    else
		      {*dptr = 0.0;}
		    dptr++;
		  }
	      }
	  }
      }

    if (k!=nnz)
      {error_msg_fatal("didn't place nnz entries in matrix");}

    /* free temp storage */
    bss_free(ptrs);
    bss_free(cols);
    bss_free(vals);
   }

/* the file is too short or number of line entries bad */
else
  {rt=ERROR; printf("\nCan't Position Someone");}

if (fclose(ifp))
  {error_msg_fatal("Can't close input file!!!");}

id->gs_handle = gs;
id->lda = lda;
id->sda = sda;
id->nnz = nnz;
id->permutation = NULL;
id->row_lens = NULL;
id->row_indices = NULL;
id->row_vals = NULL;
/*
id->col_lens = NULL;
id->col_indices = NULL;
id->col_vals = NULL;
*/

return(rt);
}




/********************************sparse_matrix.c*******************************
Function: SMI_rsb

Input : 
Output: 
Return: 
Description:  
*********************************sparse_matrix.c******************************/
int
SMI_csc(int *ptrs, int *cols, double *vals, int lda, int type)
{
  int i, j, k, rt;
  matrix_id_ptr id, SMI_get_matrix(int *);
  

  /* new template and add to list */
  if ((id = SMI_get_matrix(&rt)) == NULL)
    {error_msg_fatal("No matrix handles left!!! Increase MAX_MAT.");}

  /* initialize adj list */
  id->adj_list = init_adj_list(lda);

  /* values? */
  if (!vals)
    {error_msg_warning("SMI_csc() :: No a(i,j)'s.\n");}

  /* place i,j pairs into adj list */
  /* symmetric case but don't add diags twice */
  if (type==SYMM)
    {
      if (vals)
	{
	  for (k=i=0; i<lda; i++)
	    {
	      for (j=0; j<(*(ptrs+i+1)-*(ptrs+i)); j++)
		{
		  insert_adj_list(id->adj_list,i,*(cols+k)-1,*(vals+k)); 
		  if (i!=*(cols+k)-1)
		    {insert_adj_list(id->adj_list,*(cols+k)-1,i,*(vals+k));}
		  k++;
		}
	    }
	}
      else
	{
	  for (k=i=0; i<lda; i++)
	    {
	      for (j=0; j<(*(ptrs+i+1)-*(ptrs+i)); j++)
		{
		  insert_adj_list(id->adj_list,i,*(cols+k)-1,0.0); 
		  if (i!=*(cols+k)-1)
		    {insert_adj_list(id->adj_list,*(cols+k)-1,i,0.0);}
		  k++;
		}
	    }
	}
    }
  /* type == FULL */
  else
    {
      if (vals)
	{
	  for (k=i=0; i<lda; i++)
	    {
	      for (j=0; j<(*(ptrs+i+1)-*(ptrs+i)); j++)
		{
		  insert_adj_list(id->adj_list,i,*(cols+k)-1,*(vals+k)); 
		  k++;
		}
	    }
	}
      else
	{
	  for (k=i=0; i<lda; i++)
	    {
	      for (j=0; j<(*(ptrs+i+1)-*(ptrs+i)); j++)
		{
		  insert_adj_list(id->adj_list,i,*(cols+k)-1,0.0); 
		  k++;
		}
	    }
	}
    }

  id->lda = lda;
  id->sda = lda;
  id->nnz = k;
  id->permutation = NULL;
  id->row_vals = NULL;
  id->col_vals = NULL;

  return(rt);
}




/********************************sparse_matrix.c*******************************
Function: SMI_rsb

Input : 
Output: 
Return: 
Description:  
*********************************sparse_matrix.c******************************/
int
SMI_csc_map(int *map, int *ptrs, int *cols, double *vals, int lda, int type)
{
  int i, j, k, rt, i_map, j_map;
  matrix_id_ptr id, SMI_get_matrix(int *);
  int *inv_map, *tmp_map;


  /* new template and add to list */
  if ((id = SMI_get_matrix(&rt)) == NULL)
    {error_msg_fatal("No matrix handles left!!! Increase MAX_MAT.");}

  /* initialize adj list */
  id->adj_list = init_adj_list(lda);

  /* values? */
#ifdef INFO
  if (!vals)
    {printf("SMI_csc() :: No a(i,j)'s.");}
#endif

  tmp_map = (int *) bss_malloc(lda*INT_LEN);
  ivec_copy(tmp_map,map,lda);
  inv_map = (int *) bss_malloc(lda*INT_LEN);
  ivec_c_index(inv_map,lda);
  SMI_sort(tmp_map, inv_map, lda, SORT_INTEGER);

  /*
  printf("\n");
  for (i=0;i<lda;i++)
    {printf("%2d ",inv_map[i]);}
  printf("\n");
  */

  /* place i,j pairs into adj list */
  /* symmetric case but don't add diags twice */
  if (type==SYMM)
    {
      if (vals)
	{
	  for (k=i=0; i<lda; i++)
	    {
	      for (j=0; j<(*(ptrs+i+1)-*(ptrs+i)); j++)
		{
		  insert_adj_list(id->adj_list,i,*(cols+k)-1,*(vals+k)); 
		  if (i!=*(cols+k)-1)
		    {insert_adj_list(id->adj_list,*(cols+k)-1,i,*(vals+k));}
		  k++;
		}
	    }
	}
      else
	{
	  for (k=i=0; i<lda; i++)
	    {
	      for (j=0; j<(*(ptrs+i+1)-*(ptrs+i)); j++)
		{
		  insert_adj_list(id->adj_list,i,*(cols+k)-1,0.0); 
		  if (i!=*(cols+k)-1)
		    {insert_adj_list(id->adj_list,*(cols+k)-1,i,0.0);}
		  k++;
		}
	    }
	}
    }
  /* type == FULL */
  else
    {
      if (vals)
	{
	  /* printf("\n"); */
	  for (k=i=0; i<lda; i++)
	    {
	      i_map=inv_map[i];
	      for (j=0; j<(*(ptrs+i+1)-*(ptrs+i)); j++)
		{
		  j_map=inv_map[*(cols+k)-1];
		  insert_adj_list(id->adj_list,i_map,j_map,*(vals+k));
		  /* printf("(%d,%d)->(%d,%d) ",i+1,*(cols+k),i_map,j_map); */
		  k++;
		}
	      /* printf("\n"); */
	    }
	}
      else
	{
	  /* printf("\n"); */
	  for (k=i=0; i<lda; i++)
	    {
	      i_map=inv_map[i];
	      for (j=0; j<(*(ptrs+i+1)-*(ptrs+i)); j++)
		{
		  j_map=inv_map[*(cols+k)-1];
		  insert_adj_list(id->adj_list,i_map,j_map,0.0);
		  /* printf("(%d,%d)->(%d,%d) ",i+1,*(cols+k),i_map,j_map); */
		  k++;
		}
	      /* printf("\n"); */
	    }
	}
    }

  id->lda = lda;
  id->sda = lda;
  id->nnz = k;
  id->permutation = NULL;
  id->row_vals = NULL;
  id->col_vals = NULL;
  bss_free(tmp_map);
  bss_free(inv_map);
  
  return(rt);
}

/*************************************************************************** 
Function gaussian_elim(): solves the system Az=b for z using gaussian 
elimination modified for tri-diagonal matrices which runs in O(n) time. 
All other components of the program run in O(n) time -> the entire program
runs in O(n) time. Recall n is the number of intervals in the partition. 
This function assumes that b is passes in through the vector z.
***************************************************************************/
void
tri_calc_evec(double *d, double *sd, int n, double lambda, double *vec)
{
  int i, off;
  double *tri;
  double beta;


  /* allocate space for eigenvecto and tridiag matrix */
  dvec_zero(vec,(n+1));
  tri = (double *) bss_malloc(4*n*sizeof(double));
  dvec_zero(tri,3*n);

  /* create tri from diag and sub-diag */
  off = 0;
  *(tri+off++) = 0.0;
  *(tri+off++) = *(d+1) - lambda;
  *(tri+off)   = *(sd+2);
  for(i=1; i<(n-1); i++)
    {
      off = 3*i;
      *(tri+off++) = *(sd+i+1);
      *(tri+off++) = *(d+i+1) - lambda;
      *(tri+off)   = *(sd+i+2);
    }
  off = 3*(n-1);
  *(tri+off++) = *(sd+n);
  *(tri+off++) = *(d+n) - lambda;
  *(tri+off)   = 0.0;

  /* Down */
  for(i=0;i<(n-1);i++)
    {
      off = 3*i;
      /*
      printf("%d  %.16e\n",i,*(tri+off+1));
      printf("%d  %.16e\n",i,*(tri+off+4));
      */
      beta = *(tri+off+1);
      if (fabs(beta)<1.0e-10)
	{printf("%d  %.16e\n",i,beta);}
      *(tri+off+4) -= *(tri+off+2) * *(tri+off+3)/beta;
      /*
      if (i==(n-2))
	{printf("%d  %.16e\n",n,*(tri+off+4));}
	*/
    }

  i=(n-2);
  if (fabs(*(tri+off+4))>1.0e-7)
    {
      vec[0]=0.0;
      dvec_set(vec+1,1.0,n+1);
      printf("tri_calc_evec failure! %d ::  %.16e\n",n,*(tri+off+4));
      normalize(vec+1,n);
      return;
    }

  for(i=0;i<(n-1);i++)
    {
      off = 3*i;
      *(tri+off+2) /= *(tri+off+1);
    }

  vec[n-2]=*(tri+3*(n-2)+2);
  vec[n-1]=-1.0;

  /* up */
  for(i=(n-2);i>0;i--)
    {*(vec+i-1) = -*(tri+3*i-1) * *(vec+i);}

  normalize(vec,n);
  for(i=n;i>0;i--)
    {vec[i]=-1.0*vec[i-1];}
  vec[0]=0.0;

  bss_free(tri);
}


/*************************************************************************** 
Function gaussian_elim(): solves the system Az=b for z using gaussian 
elimination modified for tri-diagonal matrices which runs in O(n) time. 
All other components of the program run in O(n) time -> the entire program
runs in O(n) time. Recall n is the number of intervals in the partition. 
This function assumes that b is passes in through the vector z.
***************************************************************************/
int
trid_slv(double *d, double *sd, int n, double lambda, double *b)
{
  int i, off;
  double *vec, *tri;
  double beta;
  double div;

  for(i=0;i<n;i++)
    {b[i]=b[i+1];}
  b[n]=0.0;

  /* allocate space for eigenvector and tridiag matrix */
  tri = (double *) bss_malloc(4*n*sizeof(double));
  dvec_zero(tri,3*n);

  /* create tri from diag and sub-diag */
  off = 0;
  *(tri+off++) = 0.0;
  *(tri+off++) = *(d+1) - lambda;
  *(tri+off)   = *(sd+2);
  for(i=1; i<(n-1); i++)
    {
      off = 3*i;
      *(tri+off++) = *(sd+i+1);
      *(tri+off++) = *(d+i+1) - lambda;
      *(tri+off)   = *(sd+i+2);
    }
  off = 3*(n-1);
  *(tri+off++) = *(sd+n);
  *(tri+off++) = *(d+n) - lambda;
  *(tri+off)   = 0.0;


  /* Down */
  for(i=0;i<(n-1);i++)
    {
      div = *(tri+3*i+1);
      if (div==0.0)
	{
#ifdef DEBUG
	  printf("trid_slv() %d :: div=%.16e\n",i,div);
#endif
	  bss_free(tri);
	  return(0);
	}

     *(tri+3*i+2) /= div;
     *(b+i) /= div;
     *(b+i+1) += -*(tri+3*i+3)* *(b+i);
     *(tri+3*i+4) += -*(tri+3*i+3)* *(tri+3*i+2);
   }

  /* last entry */
  if (n>1)
    {
      div = *(tri+3*(n-1)+1);
      if (div==0.0)
	{
#ifdef DEBUG
	  printf("trid_slv() %d :: div=%.16e\n",i,div);
#endif
	  bss_free(tri);
	  return(0);
	}
      *(b+n-1) /= div;
    }

  /* up */
  for(i=(n-1);i>0;i--)
    {*(b+i-1) += -*(tri+3*i-1) * *(b+i);}


  normalize(b,n);
  for(i=n;i>0;i--)
    {b[i]=b[i-1];}
  b[0]=0.0;

  bss_free(tri);
  return(1);
}



/*************************************************************************** 
Function gaussian_elim(): solves the system Az=b for z using gaussian 
elimination modified for tri-diagonal matrices which runs in O(n) time. 
All other components of the program run in O(n) time -> the entire program
runs in O(n) time. Recall n is the number of intervals in the partition. 
This function assumes that b is passes in through the vector z.
***************************************************************************/
double *
trid_mult(double *d, double *sd, int n, double *b)
{
  int i, off;
  double *tmp2, *tmp;


  tmp2 = tmp = (double *)bss_malloc((2*n+1)*sizeof(double));
  dvec_zero(tmp,n+1);
  b++;
  tmp++;
  

  tmp[0] +=  d[1]*b[0];
  tmp[0] += sd[2]*b[1];
  for(i=1;i<(n-1);i++)
    {
      tmp[i] += sd[i+1]*b[i-1];
      tmp[i] +=  d[i+1]*b[i];
      tmp[i] += sd[i+2]*b[i+1];
    }
  tmp[i] += sd[i+1]*b[i-1];
  tmp[i] +=  d[i+1]*b[i];


  return(tmp2);
}



/********************************sparse_matrix.c*******************************
Function: lanczos()

Input : 
Output: 
Return: 
Description:  
*********************************sparse_matrix.c******************************/
double *
cg(double *rhs, int n, int *ia, int *ja, double *vals, double lambda)
{
  int i,j,k=0;
  double *r_k=NULL, *p_k=NULL, *x_k=NULL, *Ap_k=NULL, *r_k_old, *r_k_new;
  double *Ar_k=NULL;
  double alpha=0.0, beta=0.0, rr=0.0, rr_limit, new_rr;
  double alpha_old=0.0, beta_old=0.0;
  int *color, *ptr;
  int iter_limit;
  queue_ADT ln_q;
  struct lanczos_node *ln_ptr;
  double  *dv, *ev, *dv_w, *ev_w, dmax, dmin, **z, **Q;
  double *w_k, *tmp, *v_k, *v_k_old;
  double la, lb, lb_old;
  int *companion;
  double pAp=0.0, pAp_old=0.0;
  double rr_old;


  /*
  hmt_ortho(rhs,n);
  normalize(rhs,n);
  */

  /* initial rhs w/2-norm of ~1.0 and orthogonal to (1,1,...,1)^T */
  r_k     = (double *) bss_malloc(n*sizeof(double));
  r_k_old = (double *) bss_malloc(n*sizeof(double));
  r_k_new = (double *) bss_malloc(n*sizeof(double));
  p_k = (double *) bss_malloc(n*sizeof(double));
  dvec_copy(r_k,rhs,n);
  /*
  hmt_ortho(r_k,n);
  normalize(r_k,n);
  */
  dvec_copy(p_k,r_k,n);
  rr  = hmt_ddot(r_k,r_k,n);

  /* x_k to initial guess (the 0.0 vector) */
  x_k = (double *) bss_malloc(n*sizeof(double));
  dvec_zero(x_k,n);
  /* dvec_copy(x_k,r_k,n); */

  /* space for Ap_k */
  Ap_k = (double *) bss_malloc(n*sizeof(double));
  dvec_zero(Ap_k,n);
  Ar_k = (double *) bss_malloc(n*sizeof(double));
  dvec_zero(Ar_k,n);
  
  /* set lanczos limits */
  /*rr_limit = sqrt(rr)*RES_LIMIT;*/
  rr_limit = res_limit;
  iter_limit = MIN(n,CG_LIMIT);


  /*while ((rr>1.0e-14)&&(k<iter_limit))*/
#ifdef INFO  
  printf("\n\n");
#endif
  i=0;
  /* while ((sqrt(rr)>rr_limit) && (k<iter_limit)) */
  while ((rr>rr_limit) && (k<iter_limit)) 
  {

    printf("cg 2-norm r_k[%d]=%.16e\n",k,rr);

    k++;
    dvec_copy(r_k_old,r_k,n);
    normalize(r_k_old,n);

    skew_matrix_mult(Ap_k,p_k,n,ia,ja,vals,lambda);

    pAp_old = pAp;
    pAp = hmt_ddot(Ap_k,p_k,n);
    printf("[%d]= (Ap_k,p_k)=%.16e\n",k,pAp);
    alpha = rr/hmt_ddot(Ap_k,p_k,n);
    /* fdaxpy(n,alpha,p_k,1,x_k,1); */
    hmt_daxpy(x_k,alpha,p_k,n);

    /* fdaxpy(n,-1.0*alpha,Ap_k,1,r_k,1); */
    hmt_daxpy(r_k,-1.0*alpha,Ap_k,n);
    hmt_ortho(r_k,n);

    dvec_copy(r_k_new,r_k,n);
    normalize(r_k_new,n);

    new_rr = hmt_ddot(r_k,r_k,n);
    beta = new_rr/rr;
    rev_daxpy(r_k,beta,p_k,n);

    rr_old = rr;
    rr = new_rr;
    printf("[%d] a=%.16e\n",k,alpha);
    printf("[%d] b=%.16e\n",k,beta);
    alpha_old = alpha;
    beta_old = beta;
  }
  printf("cg 2-norm r_k[%d]=%.16e\n",k,rr);

  bss_free(p_k);
  bss_free(r_k);
  bss_free(r_k_old);
  bss_free(r_k_new);

  return(x_k);
}


/********************************sparse_matrix.c*******************************
Function: lanczos()

Input : 
Output: 
Return: 
Description:  
*********************************sparse_matrix.c******************************/
double *
cg_min_res(double *rhs, int n, int *ia, int *ja, double *vals, double lambda)
{
  int i,j,k=0;
  double *r_k=NULL, *p_k=NULL, *x_k=NULL, *Ap_k=NULL, *r_k_old, *r_k_new;
  double *Ar_k=NULL;
  double alpha=0.0, beta=0.0, rAr=0.0, new_rAr, rr_limit, rr, new_rr;
  double alpha_old=0.0, beta_old=0.0;
  int *color, *ptr;
  int iter_limit;
  queue_ADT ln_q;
  struct lanczos_node *ln_ptr;
  double  *dv, *ev, *dv_w, *ev_w, dmax, dmin, **z, **Q;
  double *w_k, *tmp, *v_k, *v_k_old;
  double la, lb, lb_old;
  int *companion;
  double pAp=0.0, pAp_old=0.0;
  double rr_old;

  /*
  hmt_ortho(rhs,n);
  normalize(rhs,n);
  */

  /* initial rhs w/2-norm of ~1.0 and orthogonal to (1,1,...,1)^T */
  r_k     = (double *) bss_malloc(n*sizeof(double));
  r_k_old = (double *) bss_malloc(n*sizeof(double));
  r_k_new = (double *) bss_malloc(n*sizeof(double));
  p_k = (double *) bss_malloc(n*sizeof(double));
  dvec_copy(r_k,rhs,n);
  normalize(r_k,n);
  /*
  hmt_ortho(r_k,n);
  */
  dvec_copy(p_k,r_k,n);

  /* x_k to initial guess (the 0.0 vector) */
  x_k = (double *) bss_malloc(n*sizeof(double));
  dvec_zero(x_k,n);
  /*  dvec_copy(x_k,r_k,n); */

  /* space for Ap_k */
  Ap_k = (double *) bss_malloc(n*sizeof(double));
  dvec_zero(Ap_k,n);
  Ar_k = (double *) bss_malloc(n*sizeof(double));
  dvec_zero(Ar_k,n);
  
  /* set lanczos limits */
  /*rr_limit = sqrt(rr)*RES_LIMIT;*/
  /*rAr_limit = res_limit;*/
  rr_limit = 1.0e-15;
  iter_limit = MIN(n,CG_MR_LIMIT);


  /*while ((rr>1.0e-14)&&(k<iter_limit))*/
#ifdef INFO  
  printf("\n\n");
#endif
  i=0;
  skew_matrix_mult(Ar_k,r_k,n,ia,ja,vals,lambda);
  rAr = hmt_ddot(r_k,Ar_k,n);
  rr = hmt_ddot(r_k,r_k,n);
  /*
  printf("res_limit=%.16e\n",res_limit);
  printf("rAr_limit=%.16e\n",rAr_limit);
  printf("rAr 2-norm r_k[start]=%.16e\n",rAr);
  */
  /*while ((fabs(rAr)>rAr_limit) && (k<iter_limit)) */
  while ((rr>rr_limit) && (k<iter_limit))
    {
      /*printf("cg 2-norm (r_k,r_k)[%d]=%.16e\n",k,rr);*/
      skew_matrix_mult(Ap_k,p_k,n,ia,ja,vals,lambda);
      
      alpha = rAr/hmt_ddot(Ap_k,Ap_k,n);
      hmt_daxpy(x_k,alpha,p_k,n);    
      hmt_daxpy(r_k,-1.0*alpha,Ap_k,n);
      rr = hmt_ddot(r_k,r_k,n);
      skew_matrix_mult(Ar_k,r_k,n,ia,ja,vals,lambda);
      new_rAr = hmt_ddot(r_k,Ar_k,n);
      beta = new_rAr/rAr;
      rev_daxpy(r_k,beta,p_k,n);
      rAr = new_rAr;
      /*
      printf("[%d] a=%.16e\n",k,alpha);
      printf("[%d] b=%.16e\n",k,beta);
      */
      k++;
    }
  printf("cg_min_res 2-norm r_k[%d]=%.16e\n",k,rAr);

  bss_free(p_k);
  bss_free(r_k);
  bss_free(r_k_old);
  bss_free(r_k_new);

  return(x_k);
}


/********************************ivec.c**************************************
Function dvec_zero()

Input : 
Output: 
Return: 
Description: 
*********************************ivec.c*************************************/
static
void 
dvec_zero(register double *arg1, register int n)
{
  while (n--)  {*arg1++ = 0.0;}
}


/********************************ivec.c**************************************
Function dvec_copy()

Input : 
Output: 
Return: 
Description: 
*********************************ivec.c*************************************/
static
void 
dvec_copy(register double *arg1, register double *arg2, register int n)
{
  while (n--)  {*arg1++ = *arg2++;}
}


/**********************************ivec.c**************************************
Function dvec_set()

Input : 
Output: 
Return: 
Description: 
***********************************ivec.c*************************************/
static
void
dvec_set(register double *arg1, register double arg2, register int n)
{
  while (n--)  {*arg1++ = arg2;}
}


/**********************************ivec.c**************************************
Function dvec_dump

Input : 
Output: 
Return: 
Description: 
***********************************ivec.c*************************************/
void
dvec_dump(double *v, int n, int tag, int tag2, char * s)
{
  int i;
  printf("%2d %2d %s :: ",tag,tag2,s);
  for (i=0;i<n;i++)
    {printf("%f ",v[i]);}
  printf("\n");
  fflush(stdout);
}


/******************************************************************************
Function: my_sort().
Input : offset of list to be sorted, number of elements to be sorted.
Output: sorted list (in ascending order).
Return: none.
Description: stack based (nonrecursive) quicksort w/brute-shell bottom. 
******************************************************************************/
static
void
dvec_sort(register double *ar, register int Size)
{
  register double *pi, *pj, temp;
  register double **top_a = (double **)doffset_stack;
  register PTRINT *top_s = pdsize_stack, *bottom_s = pdsize_stack; 
  register PTRINT size = (PTRINT) Size;

  /* we're really interested in the offset of the last element */
  /* ==> length of the list is now size + 1                    */
  size--;

  /* do until we're done ... return when stack is exhausted */
  for (;;)
    {
      /* if list is large enough use quicksort partition exchange code */
      if (size > DSORT_OPT)
	{	
	  /* start up pointer at element 1 and down at size     */  
	  pi = ar+1;
	  pj = ar+size;

	  /* find middle element in list and swap w/ element 1 */
	  SWAP(*(ar+(size>>1)),*pi)

	  pj = ar+size; 

	  /* order element 0,1,size-1 st {M,L,...,U} w/L<=M<=U */
	  /* note ==> pivot_value in index 0                   */
	  if (*pi > *pj)  
	    {SWAP(*pi,*pj)}
	  if (*ar > *pj) 
	    {SWAP(*ar,*pj)}
	  else if (*pi > *ar)   
	    {SWAP(*(ar),*(ar+1))}

	  /* partition about pivot_value ...  	                    */
	  /* note lists of length 2 are not guaranteed to be sorted */
	  for(;;)
	    {
	      /* walk up ... and down ... swap if equal to pivot! */
	      do pi++; while (*pi<*ar);
	      do pj--; while (*pj>*ar);

	      /* if we've crossed we're done */
	      if (pj<pi) break;

	      /* else swap */
	      SWAP(*pi,*pj)
	    }

	  /* place pivot_value in it's correct location */
	  SWAP(*ar,*pj)

	  /* test stack_size to see if we've exhausted our stack */
	  if (top_s-bottom_s >= DSORT_STACK)
	    {error_msg_fatal("\nSTACK EXHAUSTED!!!\n");}

	  /* push right hand child iff length > 1 */
	  if ((*top_s = size-(pi-ar)))
	    {
	      *(top_a++) = pi;
	      size -= *top_s+2;  
	      top_s++;
	    }
	  /* set up for next loop iff there is something to do */
	  else if (size -= *top_s+2) 
	    {;}
	  /* might as well pop - note NR_OPT >=2 ==> we're ok! */
	  else
	    {
	      ar = *(--top_a);
	      size = *(--top_s);
	    }
	}

      /* else sort small list directly then pop another off stack */
      else
	{
	  /* insertion sort for bottom */
          for (pj=ar+1;pj<=ar+size;pj++)
            {
              temp = *pj;
              for (pi=pj-1;pi>=ar;pi--)
                {
                  if (*pi <= temp) break;
                  *(pi+1)=*pi;
                }
              *(pi+1)=temp;
	    }

	  /* check to see if stack is exhausted ==> DONE */
	  if (top_s==bottom_s) return;
	  
	  /* else pop another list from the stack */
	  ar = *(--top_a);
	  size = *(--top_s);
	}
    }
}


/******************************************************************************
Function: my_sort().
Input : offset of list to be sorted, number of elements to be sorted.
Output: sorted list (in ascending order).
Return: none.
Description: stack based (nonrecursive) quicksort w/brute-shell bottom. 
******************************************************************************/
static
void
dvec_sort_companion(register double *ar, register int *ar2, register int Size)
{
  register double *pi, *pj, temp;
  register double **top_a = (double **)doffset_stack;
  register PTRINT *top_s = pdsize_stack, *bottom_s = pdsize_stack; 
  register PTRINT size = (PTRINT) Size;

  register int *pi2, *pj2;
  register int ptr;
  register PTRINT mid;


  /* we're really interested in the offset of the last element */
  /* ==> length of the list is now size + 1                    */
  size--;

  /* do until we're done ... return when stack is exhausted */
  for (;;)
    {
      /* if list is large enough use quicksort partition exchange code */
      if (size > DSORT_OPT)
	{	
	  /* start up pointer at element 1 and down at size     */  
	  mid = size>>1;
	  pi = ar+1;
	  pj = ar+mid;
	  pi2 = ar2+1;
	  pj2 = ar2+mid;

	  /* find middle element in list and swap w/ element 1 */
	  SWAP(*pi,*pj)
	  P_SWAP(*pi2,*pj2)

	  /* order element 0,1,size-1 st {M,L,...,U} w/L<=M<=U */
	  /* note ==> pivot_value in index 0                   */
	  pj = ar+size;
	  pj2 = ar2+size;
	  if (*pi > *pj)  
	    {SWAP(*pi,*pj) P_SWAP(*pi2,*pj2)}
	  if (*ar > *pj) 
	    {SWAP(*ar,*pj) P_SWAP(*ar2,*pj2)}
	  else if (*pi > *ar)   
	    {SWAP(*(ar),*(ar+1)) P_SWAP(*(ar2),*(ar2+1))}

	  /* partition about pivot_value ...  	                    */
	  /* note lists of length 2 are not guaranteed to be sorted */
	  for(;;)
	    {
	      /* walk up ... and down ... swap if equal to pivot! */
	      do {pi++; pi2++;} while (*pi<*ar);
	      do {pj--; pj2--;} while (*pj>*ar);

	      /* if we've crossed we're done */
	      if (pj<pi) break;

	      /* else swap */
	      SWAP(*pi,*pj)
	      P_SWAP(*pi2,*pj2)
	    }

	  /* place pivot_value in it's correct location */
	  SWAP(*ar,*pj)
	  P_SWAP(*ar2,*pj2)

	  /* test stack_size to see if we've exhausted our stack */
	  if (top_s-bottom_s >= DSORT_STACK)
	    {error_msg_fatal("\nSTACK EXHAUSTED!!!\n");}

	  /* push right hand child iff length > 1 */
	  if ((*top_s = size-(pi-ar)))
	    {
	      *(top_a++) = pi;
	      *(top_a++) = (double *) pi2;
	      size -= *top_s+2;  
	      top_s++;
	    }
	  /* set up for next loop iff there is something to do */
	  else if (size -= *top_s+2) 
	    {;}
	  /* might as well pop - note NR_OPT >=2 ==> we're ok! */
	  else
	    {
	      ar2 = (int *) *(--top_a);
	      ar  = *(--top_a);
	      size = *(--top_s);
	    }
	}

      /* else sort small list directly then pop another off stack */
      else
	{
	  /* insertion sort for bottom */
          for (pj=ar+1, pj2=ar2+1;pj<=ar+size;pj++,pj2++)
            {
              temp = *pj;
              ptr = *pj2;
              for (pi=pj-1,pi2=pj2-1;pi>=ar;pi--,pi2--)
                {
                  if (*pi <= temp) break;
                  *(pi+1)=*pi;
                  *(pi2+1)=*pi2;
                }
              *(pi+1)=temp;
              *(pi2+1)=ptr;
	    }

	  /* check to see if stack is exhausted ==> DONE */
	  if (top_s==bottom_s) return;
	  
	  /* else pop another list from the stack */
	  ar2 = (int *) *(--top_a);
	  ar  = *(--top_a);
	  size = *(--top_s);
	}
    }
}
