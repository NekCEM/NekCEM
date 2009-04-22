/** \file tools.c
    \brief Tools related to XXT, such as moving window element to processor mapping 
    
   	A detailed description goes here 
*/

/*************************************xxt.c************************************
Module Name: xxt
Module Info:

Author:  Henry M. Tufo III

e-mail:  hmt@cs.brown.edu

sn-mail: Division of Applied Mathematics,
	 Brown University,
         Box F
	 Providence, RI 02912

Tel:	 (401) 863-7666


Last Modification: 10.26.97
**************************************xxt.c***********************************/


/*************************************xxt.c************************************
NOTES ON USAGE:

**************************************xxt.c***********************************/


/*************************************xxt.c************************************
FILE FORMAT:
------------------------------ Begin File -------------------------------------

------------------------------ End   File -------------------------------------

Note:
**************************************xxt.c***********************************/

/* C modules for I/O etc. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include <unistd.h>

#if   defined NXSRC
#ifndef DELTA
#include <nx.h>
#endif

#ifdef DELTA
#include <mesh.h>
#include <cveclib.h>
#endif

#elif defined MPISRC
#include <mpi.h>

#endif

/* mine : const before types! */
#include "const.h"
#include "types.h"
#include "comm.h"
#include "error.h"
#include "ivec.h"
#include "bss_malloc.h"
#include "queue.h"


#define MAX_REA_NAME   120
#define STD_READ_BUF   1000

static void set_file_names(void);

/* rea file handle */
char dir_name[MAX_REA_NAME+5];
char rea_name[MAX_REA_NAME+5];
char map_name[MAX_REA_NAME+5];
char sep_name[MAX_REA_NAME+5];

/* number of spectral elements (SE) over all P and # I own */
static int nel_global, nel_local;

/* number of dof w/n_global = nvu-nvo */
/* number of unique vertices I own = n - # outflow I own */
/* after separators read in equal to E coarse dimension n */
static int n_global, n_local;

/* depth of separator tree and max number of processors we can run on */
static int depth, max_proc;

/* number of corners a SE has ... 4 or 8 only! */
static int nc;

/* number of vertices = nel_global*nc */
static int nv;

/* number of unique vertices = nv - {copies over all SE} */
static int nvu;

/* number of outflow vertices */
static int nvo;

/* holds element vertex global numbering - size=n_global*nv */
static int *vertex=NULL;

/* holds element to processor map */
static int *map;

int slices[256];	/* TODO: this is for debugging only so remove it when I get done */

static int numlines;

/*************************************xxt.c************************************
Function:

Input :
Output:
Return:
Description:
**************************************xxt.c***********************************/
#ifdef UPCASE
void
FXXT_IVERTEX_MAP  (int *x, int *nelv, int *ncr)
#elif  IBM
void
fxxt_ivertex_map  (int *x, int *nelv, int *ncr)
#else
void
fxxt_ivertex_map_ (int *x, int *nelv, int *ncr)
#endif
{
    int i,j,k;
    int *iptr_m, *iptr_v;;


#ifdef DEBUG
    error_msg_warning("xxt_ivertex_map() :: begin\n");
#endif

    if (*nelv != nel_global) {
        error_msg_fatal("nelv=%d != nel_global=%d!\n",*nelv,nel_global); return;
    }

    if (*ncr != nc) {
        error_msg_fatal("nc=%d, *ncr=%d\n",nc,*ncr); return;
    }

    /*printf("writing x:\n");*/
    /* ok ... process it */
    for (iptr_v=vertex, iptr_m=map, i=0; i<nel_global; i++) {
        for (j=0; j<nc; j++)
            *x++ = *iptr_v++;
    }

#ifdef DEBUG
    error_msg_warning("xxt_ivertex_map() :: end\n");
#endif
}

#ifdef UPCASE
void
FXXT_IVERTEX_MAPW  (int *x, int *nelv, int *ncr)
#elif  IBM
void
fxxt_ivertex_mapw  (int *x, int *nelv, int *ncr)
#else
void
fxxt_ivertex_mapw_ (int *x, int *nelv, int *ncr)
#endif
{
    int i,j,k;
    int *iptr_m, *iptr_v;;


#ifdef DEBUG
    error_msg_warning("xxt_ivertex_map() :: begin\n");
#endif

    if (*nelv != nel_global) {
        error_msg_fatal("nelv=%d != nel_global=%d!\n",*nelv,nel_global); return;
    }

    if (*ncr != nc) {
        error_msg_fatal("nc=%d, *ncr=%d\n",nc,*ncr); return;
    }

    printf("writing x:\n");
    /* ok ... process it */
    for (iptr_v=vertex, i=0; i<numlines; i++) {
        if (my_id==0) printf("%dth line: ",i);
        for (j=0; j<nc; j++) {
            *x++ = *iptr_v++;
              if (my_id==0)  printf("%d ",*(x-1));
        }
         if (my_id==0) printf("\n"); 
    }

#ifdef DEBUG
    error_msg_warning("xxt_ivertex_map() :: end\n");
#endif
}




/** Element to processor map (non moving window)
 *	
 *	Read in element to processor map information for non-moving window case. In this case, we use the mapping provided by genmap without modification.
 *
 *	@param out_map Element mapping. Position-based. I.e. out_map[0]=5 means element 0 is assigned to rank 5
 *	@param nelgt number of global elements
 *	@param dim number of vertices per element
 *
 */
#ifdef UPCASE
void
XXT_ELM_TO_PROC  (int *out_map, int *nelgt, int *dim)
#elif  IBM
void
xxt_elm_to_proc  (int *out_map, int *nelgt, int *dim)
#else
void
xxt_elm_to_proc_ (int *out_map, int *nelgt, int *dim)
#endif
{
    int  i,j,k=0;
    int  who;
    char *buf, *token;
    int  *iptr_m, *iptr_v;
    FILE *ifp;
    int  proc_number;

#ifdef DEBUG
    error_msg_warning("xxt_elm_to_proc() :: begin\n");
#endif

    set_file_names();

    if ((ifp=fopen(map_name,"r"))==NULL) {
        error_msg_fatal("can't open %s for input!\n",map_name); return;
    }

    /* place to store information from xxt elm2proc file */
    buf = (char *) bss_malloc(STD_READ_BUF*sizeof(char));

    /* read header line */
    if (fgets(buf,STD_READ_BUF,ifp) == NULL) {
        error_msg_fatal("where's the damn header line!\n");
    }

    if ((token=(char *) strtok(buf,DELIM)) == NULL) {
        error_msg_fatal("nel_global data missing!\n");
    }
    nel_global = atoi(token);

    if (nel_global!=*nelgt) {
        error_msg_fatal("nel_global=%d != nelgt=%d!\n",nel_global,*nelgt);
    }

    if ((token=(char *) strtok(NULL,DELIM)) == NULL) {
        error_msg_fatal("n_global (dof) data missing!\n");
    }
    n_global = atoi(token);

    if ((token=(char *) strtok(NULL,DELIM)) == NULL) {
        error_msg_fatal("depth data missing!\n");
    }
    depth = atoi(token);

    if ((token=(char *) strtok(NULL,DELIM)) == NULL) {
        error_msg_fatal("max proc data missing!\n");
    }
    max_proc = atoi(token);

    if (num_nodes>max_proc) {
        error_msg_fatal("max_proc=%d < num_nodes=%d!\n",max_proc,num_nodes);
    }

    if ((token=(char *) strtok(NULL,DELIM)) == NULL) {
        error_msg_fatal("nv data missing!\n");
    }
    nv = atoi(token);

    if ((token=(char *) strtok(NULL,DELIM)) == NULL) {
        error_msg_fatal("nvu data missing!\n");
    }
    nvu = atoi(token);

    if ((token=(char *) strtok(NULL,DELIM)) == NULL) {
        error_msg_fatal("nvo data missing!\n");
    }
    nvo = atoi(token);

    nc = nv/nel_global;
    if (nc!=*dim) {
        error_msg_fatal("nc=%d but %d passed in!\n",nc,*dim);
    }


    /* grab space for data input */
    map    = iptr_m = (int *) bss_malloc(nel_global*INT_LEN);
    vertex = iptr_v = (int *) bss_malloc(nv*INT_LEN);

    /* map from max_proc to num_nodes */
    for (i=num_nodes, k=0; i<max_proc; k++, i<<=1) {}
    if  (i!=max_proc) {
        error_msg_fatal("%d,%d ==> num_nodes != 2^k!\n",i,max_proc);
    }

    /* read in element to processor map and rsb vertex numberings */
    for (nel_local=i=0; i<nel_global; i++) {
        /* read elm i's data */
        if (fgets(buf,STD_READ_BUF,ifp) == NULL) {
            error_msg_fatal("where's the damn %d'th line!\n",i);
        }

        /*      if (9<=i && i<18)
         *            {
         *                   printf("======%d\n",i);*/

        /* first field hold processor number in 0,...,max_proc-1 */
        if ((token=(char *) strtok(buf,DELIM)) == NULL) {
            error_msg_fatal("proc data missing!\n");
        }
        proc_number = atoi(token);
        *iptr_m++ = *out_map++ = who = atoi(token)>>k;
        /* printf("---> %d :: %d :: %d: %d\n",who,*iptr_m,*(iptr_m-1),proc_number);*/
        /* do I own it? */
        if (my_id == who) {
            nel_local++;
        }

#ifdef DEBUG
        printf("%d :: %d, %d\n",i,who,my_id);
#endif

        /* uuuuuu */
        /* remaining fields hold vertex global numbers in HC order! */
        for (j=0; j<nc; j++) {
            if ((token=(char *) strtok(NULL,DELIM)) == NULL) {
                error_msg_fatal("elm#%d vertex#%d data missing!\n",i+1,j+1);
            }
            *iptr_v++ = atoi(token);
            /*      printf(" %d ",*(iptr_v-1)); */
        }

    }

    fclose(ifp);
    bss_free(buf);

#ifdef DEBUG
    error_msg_warning("xxt_elm_to_proc() :: nel_global=%d, nel_local=%d\n",
                      nel_global,nel_local);
    error_msg_warning("xxt_elm_to_proc() :: end\n");
#endif
}


/**
 *	Read in element to processor map information for moving window case.
 *
 *	@param out_map Element mapping. Position-based. I.e. out_map[0]=5 means element 0 is assigned to rank 5
 *	@param nelgt number of global elements
 *	@param dim number of vertices per element
 *	@param start element number for start of window - 1
 *	@param end element number for end of window
 *	@param se number of elements in a slice
 *
 */
#ifdef UPCASE
void
XXT_ELM_TO_PROCW  (int *out_map, int *nelgt, int *dim, int *start, int *end, int *se)
#elif  IBM
void
xxt_elm_to_procw  (int *out_map, int *nelgt, int *dim, int *start, int *end, int *se)
#else
void
xxt_elm_to_procw_ (int *out_map, int *nelgt, int *dim, int *start, int *end, int* se)
#endif
{

    int  i,fl, j,k=0;
    char *buf, *token;
    int  *iptr_m, *iptr_v;
    FILE *ifp;
    static int  nap=0; 	           /* number of available processors */
    static int  solw=0;        	   
    static int  eolw=0;        	   
    static int  slice=0;           /* the current slice */
    int  eps = *se;		   /* elements per slice */
    int  sp, ep;	           /* starting process, ending process */
    int  ns;/* = (*end-*start+1)/ eps; number of slices in this window */
    int  pps;                      /* processors per slice */
    int  epp;		           /* elements per processor */
    int  resp=0;	           /* who is responsible for a particular element */
    int  ats=0;                    /* assigned this slice */
    int  atp=0;		           /* how many elements have we assigned to the last processor? */
	static int firstTime=1;		/* we have to do things slightly different the first time through */

    *start = *start - 1;

/*  */
    if(my_id==0)printf("ARGUMENTS: nelgt=%d dim=%d start=%d newstart=%d end=%d se=%d\n", *nelgt, *dim, *start, *end - *se, *end, *se);
/* */

/*
    printf("ARG: nelgt=%d dim=%d start=%d newstart=%d end=%d se=%d\n", *nelgt, *dim, *start, *end - *se, *end, *se);
*/
	if (!firstTime) *start = *end - *se;
    ns = (*end-*start+1)/ eps;
    printf("ARGUMENTS: end=%d start=%d eps=%d ns=%d\n", *end, *start, eps,ns);

#ifdef DEBUG
    error_msg_warning("xxt_elm_to_procw() :: begin\n");
#endif

    set_file_names();

    if ((ifp=fopen(map_name,"r"))==NULL)
        error_msg_fatal("can't open %s for input!\n",map_name);

    /* place to store information from xxt elm2proc file */
    buf = (char *) bss_malloc(STD_READ_BUF*sizeof(char));

    /* read header line */
    if (fgets(buf,STD_READ_BUF,ifp) == NULL) {
        error_msg_fatal("where's the damn header line!\n");
    }

    if ((token=(char *) strtok(buf,DELIM)) == NULL) {
        error_msg_fatal("nel_global data missing!\n");
    }
    nel_global = atoi(token);
    printf("ARG: nel_global=%d \n", nel_global);

    if (nel_global!=*nelgt) {
        error_msg_fatal("nel_global=%d != nelgt=%d!\n",nel_global,*nelgt);
    }

    if ((token=(char *) strtok(NULL,DELIM)) == NULL) {
        error_msg_fatal("n_global (dof) data missing!\n");
    }
    n_global = atoi(token);
    printf("ARG: n_global=%d \n", n_global);

    if ((token=(char *) strtok(NULL,DELIM)) == NULL) {
        error_msg_fatal("depth data missing!\n");
    }
    depth = atoi(token);

    if ((token=(char *) strtok(NULL,DELIM)) == NULL) {
        error_msg_fatal("max proc data missing!\n");
    }
    max_proc = atoi(token);
    printf("ARG: max_proc=%d \n", max_proc);

    if (num_nodes>max_proc) {
        error_msg_fatal("max_proc=%d < num_nodes=%d!\n",max_proc,num_nodes);
    }

    if ((token=(char *) strtok(NULL,DELIM)) == NULL) {
        error_msg_fatal("nv data missing!\n");
    }
    nv = atoi(token);
    printf("ARG: nv=%d \n", nv);

    if ((token=(char *) strtok(NULL,DELIM)) == NULL) {
        error_msg_fatal("nvu data missing!\n");
    }
    nvu = atoi(token);
    printf("ARG: nvu=%d \n", nvu);

    if ((token=(char *) strtok(NULL,DELIM)) == NULL) {
        error_msg_fatal("nvo data missing!\n");
    }
    nvo = atoi(token);
    printf("ARG: nvo=%d \n", nvo);

    nc = nv/nel_global;
    printf("ARG: nc=%d \n", nc);

    if (nc!=*dim) {
        error_msg_fatal("nc=%d but %d passed in!\n",nc,*dim);
    }

    /* grab space for data input */
    map = iptr_m = (int *) bss_malloc((*end-*start+1)*INT_LEN);
    printf("ARG: map=%d int_len=%d malloc=%d \n", map,INT_LEN,(*end-*start+1)*INT_LEN);

    /* vertex = iptr_v = (int *) bss_malloc((*end-*start+1)*nc*INT_LEN);*/

    if (vertex == NULL) {
        vertex = iptr_v = (int *) bss_malloc((*end-*start+1)*nc*INT_LEN);
  
    /* TODO: maybe this should use nelgt number of items, but in the future 
       the total number of elements will always be less than the number we 
       used the first time through, so this *should* be safe */
		solw=*start;			/* move to the next slice */
		eolw=*start + *se - 1;		/* mark the end as the end of the first slice. Bit of a misnomer, really */
		firstTime=1;
              printf("ARG: first time vertex=%d solw=%d eolw=%d \n", vertex,solw,eolw);
      } else {
        iptr_v = vertex; 	  /* reset the pointer back to the start */
      }    

    numlines=*end-*start+1;
    printf("ARG: numlines=%d end=%d start=%d \n", numlines,*end,*start);

    if (nap == 0) {	/* first time here, setting up the initial window */
	for(i=0;i<*nelgt;i++) out_map[i]=-2; /* clear out the map */
        nap = num_nodes;
        solw= 0;
    } else {  
        /* we've been here before... looks like we are moving the window this time */
        resp = sp = out_map[solw];     /* the first node assigned in that slice */
        ep  = out_map[eolw];    /* the last node assigned */
        nap = ep - sp + 1;		      /* set number of active processes in this case */
        printf("ARG: resp=%d solw=%d ep=%d nap=%d\n", resp,solw,ep,nap);

        /* set the old map values to -1 */
	for(i=0; i<=eolw; i++) {
        out_map[i] = -1;
        printf("ARG: i=%d eolw=%d out_map=%d\n", i,eolw,out_map[i]);
        }
    }

    pps = floor(((double)nap/ (double)ns)+0.5);                /* num of proc*/    
    epp = floor(((double)eps/ (double)nap* (double)ns) + 0.5); /*elt per proc*/
    printf("ARG: pps=%d epp=%d \n", pps,epp);
    exit(1);

/*
    if (my_id==0)printf("DEBUG: start=%d end=%d dim=%d se=%d, nelgt=%d sp=%d ep=%d nap=%d ns=%d eps=%d epp=%d pps=%d solw=%d eolw=%d\n", *start, *end, *dim, *se, *nelgt, sp, ep, nap, ns, eps, epp, pps, solw, eolw);
*/

    /* map from max_proc to num_nodes */
    
    for (i=num_nodes, k=0; i<max_proc; k++, i<<=1) {}

    if  (i!=max_proc) {
        error_msg_fatal("%d,%d ==> num_nodes != 2^k!\n",i,max_proc);
    }

    /* read through the map file */

    for (fl=0; fl < *end; fl++) {
     
        /* grab the next input line */
        if (fgets(buf,STD_READ_BUF,ifp) == NULL) {
            error_msg_fatal("reached the end of the map file unexpectedly at line %d !!\n", fl+1);
        }

        /* move through the file till we get to the window we are interested in */
        if (fl<*start) continue;  

        /* first field hold processor number in 0,...,max_proc-1 */
        /* note that we don't care about this value. The genmap code doesn't do 
           an element->processor map in a way that helps us, and at present 
           it doesn't include a way to change the default behavior */

        if ((token=(char *) strtok(buf,DELIM)) == NULL) {
            error_msg_fatal("proc data missing from line %d!\n", fl+1);
        }

        /* check if it's time for the next slice */
        if (ats == eps) {
            slice++;	        /* move to new slice */
            slices[++resp] = 1;/* used only for debugging */
            /* reset assigned counts */
            atp =0;
            ats =0;

            /* now figure out if responsibilities change. We might be a bit 
               short-staffed on the last element if it doesn't break up nice 
               and evenly, so consider that case as well */
            if (nap - resp < pps) {
                pps = nap - resp;
                epp = floor(((double)eps / (double)pps) + 0.5);
            }
        }

        /* take responsibility for the element */
        ats++;
        atp++;
        /* have we assigned all this processor is supposed to do? */
        if (atp >= epp) {
            atp=0;
            if (ats != eps){ 
				resp++;
			}
        }

		if (resp >= num_nodes){
                    if (num_nodes==1){
                        }else{
			error_msg_fatal("ERROR: walked past the end of allocation resp=%d >= num_nodes=%d\n", resp, num_nodes);
		        }
                }

        out_map[fl] = resp;
        out_map[fl] = num_nodes == 1 ? 0: resp;

        /* grab the remaining fields, which hold global vertex numbers (HC order) */
        for (j=0; j<nc; j++) {
            if ((token=(char *) strtok(NULL,DELIM)) == NULL) {
                error_msg_fatal("elm#%d vertex#%d data missing!\n",fl+1,j+1);
            }
            *iptr_v++ = atoi(token);
        }
    }

	if (firstTime){
		firstTime = 0;
	}else{
		solw += *se;			/* move to the next slice */
		eolw += *se;		/* mark the ends as the end of the first slice. Bit of a misnomer, really */
	}
    fclose(ifp);
    bss_free(buf);


#ifdef DEBUG
    error_msg_warning("xxt_elm_to_procw() :: nel_global=%d, nel_local=%d\n",
                     nel_global,nel_local);
    error_msg_warning("xxt_elm_to_procw() :: end\n");
#endif

}


/*************************************xxt.c************************************
Function:

Input :
Output:
Return:
Description:
**************************************xxt.c***********************************/
void
#ifdef UPCASE
HMT_SET_FILE_NAMES_(int *nn, char *path)
#elif  IBM
hmt_set_file_names (int *nn, char *path)
#else
hmt_set_file_names_(int *nn, char *path)
#endif
{
    int i=0;
    int pos=0;


    for (i=0; i<*nn; i++) {
#ifdef DEBUG
        printf("%c",path[i]);
#endif
        if (path[i]==' ')  {
            break;
        }
        if (path[i]=='.')  {
            pos=i;
        }
        map_name[i] = path[i];
        sep_name[i] = path[i];
        rea_name[i] = path[i];
    }

    for (i=pos+1; i<pos+4; i++) {
        map_name[pos+1] = 'm';
        map_name[pos+2] = 'a';
        map_name[pos+3] = 'p';
        map_name[pos+4] = '\0';
    }
#ifdef DEBUG
    printf("%s\n",map_name);
#endif

    for (i=pos+1; i<pos+4; i++) {
        sep_name[pos+1] = 's';
        sep_name[pos+2] = 'e';
        sep_name[pos+3] = 'p';
        sep_name[pos+4] = '\0';
    }
#ifdef DEBUG
    printf("%s\n",sep_name);
#endif

    for (i=pos+1; i<pos+4; i++) {
        rea_name[pos+4] = '\0';
    }
#ifdef DEBUG
    printf("%s\n",rea_name);
#endif

}



/*************************************xxt.c************************************
Function:

Input :
Output:
Return:
Description:
**************************************xxt.c***********************************/
void
#ifdef UPCASE
HMT_FIX_PATH(int *nn, char *path)
#elif  IBM
hmt_fix_path (int *nn, char *path)
#else
hmt_fix_path_(int *nn, char *path)
#endif
{
    int i=0;
    int pos=0;

#ifdef SINGLE_DIR
    return;
#endif

#ifdef DEBUG
    printf("%c",path[0]);
    printf("%c",path[1]);
    printf("%c",path[2]);
    printf("%c",path[3]);
    printf("%c\n",path[4]);
    printf("%d\n",*nn);
    fflush(stdout);
#endif

    for (i=0; i<*nn; i++) {
#ifdef DEBUG
        printf("%c",path[i]);
#endif
        if (path[i]==' ')  {
            break;
        }
        if (path[i]=='_')  {
            pos=i;
        }
        dir_name[i] = path[i];
    }

#ifdef DEBUG
    printf("\n");
    printf("%d\n",i);
    printf("%c\n",path[pos]);
    printf("%d\n",pos);
    fflush(stdout);
#endif

    i = (my_id+9)%16; i++;
    switch (i) {
    case 1:
        dir_name[pos+1]=path[pos+1]='1';
        dir_name[pos+2]=path[pos+2]='/';
        dir_name[pos+3]=path[pos+3]='/';
        dir_name[pos+4]='\0';
        break;
    case 2:
        dir_name[pos+1]=path[pos+1]='2';
        dir_name[pos+2]=path[pos+2]='/';
        dir_name[pos+3]=path[pos+3]='/';
        dir_name[pos+4]='\0';
        break;
    case 3:
        dir_name[pos+1]=path[pos+1]='3';
        dir_name[pos+2]=path[pos+2]='/';
        dir_name[pos+3]=path[pos+3]='/';
        dir_name[pos+4]='\0';
        break;
    case 4:
        dir_name[pos+1]=path[pos+1]='4';
        dir_name[pos+2]=path[pos+2]='/';
        dir_name[pos+3]=path[pos+3]='/';
        dir_name[pos+4]='\0';
        break;
    case 5:
        dir_name[pos+1]=path[pos+1]='5';
        dir_name[pos+2]=path[pos+2]='/';
        dir_name[pos+3]=path[pos+3]='/';
        dir_name[pos+4]='\0';
        break;
    case 6:
        dir_name[pos+1]=path[pos+1]='6';
        dir_name[pos+2]=path[pos+2]='/';
        dir_name[pos+3]=path[pos+3]='/';
        dir_name[pos+4]='\0';
        break;
    case 7:
        dir_name[pos+1]=path[pos+1]='7';
        dir_name[pos+2]=path[pos+2]='/';
        dir_name[pos+3]=path[pos+3]='/';
        dir_name[pos+4]='\0';
        break;
    case 8:
        dir_name[pos+1]=path[pos+1]='8';
        dir_name[pos+2]=path[pos+2]='/';
        dir_name[pos+3]=path[pos+3]='/';
        dir_name[pos+4]='\0';
        break;
    case 9:
        dir_name[pos+1]=path[pos+1]='9';
        dir_name[pos+2]=path[pos+2]='/';
        dir_name[pos+3]=path[pos+3]='/';
        dir_name[pos+4]='\0';
        break;
    case 10:
        dir_name[pos+1]=path[pos+1]='1';
        dir_name[pos+2]=path[pos+2]='0';
        dir_name[pos+3]=path[pos+3]='/';
        dir_name[pos+4]='\0';
        break;
    case 11:
        dir_name[pos+1]=path[pos+1]='1';
        dir_name[pos+2]=path[pos+2]='1';
        dir_name[pos+3]=path[pos+3]='/';
        dir_name[pos+4]='\0';
        break;
    case 12:
        dir_name[pos+1]=path[pos+1]='1';
        dir_name[pos+2]=path[pos+2]='2';
        dir_name[pos+3]=path[pos+3]='/';
        dir_name[pos+4]='\0';
        break;
    case 13:
        dir_name[pos+1]=path[pos+1]='1';
        /* dir_name[pos+2]=path[pos+2]='3'; */
        dir_name[pos+2]=path[pos+2]='7';
        dir_name[pos+3]=path[pos+3]='/';
        dir_name[pos+4]='\0';
        break;
    case 14:
        dir_name[pos+1]=path[pos+1]='1';
        dir_name[pos+2]=path[pos+2]='4';
        dir_name[pos+3]=path[pos+3]='/';
        dir_name[pos+4]='\0';
        break;
    case 15:
        dir_name[pos+1]=path[pos+1]='1';
        dir_name[pos+2]=path[pos+2]='5';
        dir_name[pos+3]=path[pos+3]='/';
        dir_name[pos+4]='\0';
        break;
    case 16:
        dir_name[pos+1]=path[pos+1]='1';
        dir_name[pos+2]=path[pos+2]='6';
        dir_name[pos+3]=path[pos+3]='/';
        dir_name[pos+4]='\0';
        break;
    default:
        error_msg_fatal("Ooops ... %d too large",i);
        break;
    }

#ifdef DEBUG
    printf("%s\n",dir_name);
#endif
}



/*************************************xxt.c************************************
Function:

Input :
Output:
Return:
Description:
**************************************xxt.c***********************************/
static
void
set_file_names(void) {
    char *token;
    FILE *ifp;


#ifdef DEBUG
    error_msg_warning("set_file_names() :: begin\n");
#endif

#ifdef OLD
    if ((ifp=fopen("/cacr/home/user/hmt/schwarz/nx2d10/xxt_map.rea","r"))==NULL)
        if ((ifp=fopen("/u/hmt/schwarz/nx2d10/xxt_map.rea","r"))==NULL)
            if ((ifp=fopen("xxt_map.rea","r"))==NULL)
#endif

                if ((ifp=fopen("SESSION.NAME","r"))==NULL) {
                    error_msg_fatal("can't open SESSION.NAME for input!\n"); return;
                }

    if (fgets(rea_name,MAX_REA_NAME,ifp) == NULL) {
        error_msg_fatal("SESSION.NAME empty?\n");
    }

    if (token=strpbrk(rea_name,DELIM)) {
        *token = '\0';
    }

    if (fgets(dir_name,MAX_REA_NAME,ifp) == NULL) {
        error_msg_fatal("SESSION.NAME missing line 2?\n");
    }

    if (token=strpbrk(dir_name,DELIM)) {
        *token = '\0';
    }

    strcpy(map_name,dir_name);
    strcat(map_name,rea_name);
    strcpy(sep_name,map_name);
    strcpy(rea_name,map_name);

    strcat(rea_name,".rea");
    strcat(map_name,".map");
    strcat(sep_name,".sep");

#ifdef OLD
    if (!my_id) {
        printf("%s\n",map_name);
    }


    printf("%d :: %s",strlen(rea_name),rea_name);
    printf("%d :: %s",strlen(rea_name),rea_name);
    printf("%d :: %s",strlen(map_name),map_name);
    printf("%d :: %s",strlen(sep_name),sep_name);
#endif

    fclose(ifp);

#ifdef DEBUG
    error_msg_warning("set_file_names() :: end\n");
#endif
}



/*************************************xxt.c************************************
Function:

Input :
Output:
Return:
Description:
**************************************xxt.c***********************************/
void
#ifdef UPCASE
HMT_FIX_PARAM  (REAL *htol, REAL *h2, REAL *lpc, REAL *gpc)
#elif  IBM
hmt_fix_param  (REAL *htol, REAL *h2, REAL *lpc, REAL *gpc)
#else
hmt_fix_param_ (REAL *htol, REAL *h2, REAL *lpc, REAL *gpc)
#endif

{
    char buffer[MAX_REA_NAME+1];
    char *token;
    FILE *ifp;
    double atof(const char *);


#ifdef DEBUG
    error_msg_warning("hmt_fix_parameters() :: begin\n");
#endif

    if ((ifp=fopen("SESSION.NAME","r"))==NULL) {
        error_msg_fatal("can't open SESSION.NAME for input!\n"); return;
    }

    if (fgets(buffer,MAX_REA_NAME,ifp) == NULL) {
        error_msg_fatal("SESSION.NAME empty?\n");
    }

    if (fgets(buffer,MAX_REA_NAME,ifp) == NULL) {
        error_msg_fatal("SESSION.NAME missing line 2?\n");
    }

    if (fgets(buffer,MAX_REA_NAME,ifp) == NULL) {
        error_msg_fatal("missing helmholtz tol?\n");
    }
    if (token=strpbrk(buffer,DELIM)) {
        *token = '\0';
    }
    *htol = atof(buffer);
    printf("%s,htol=%f\n",buffer,*htol);

    if (fgets(buffer,MAX_REA_NAME,ifp) == NULL) {
        error_msg_fatal("missing h2?\n");
    }
    if (token=strpbrk(buffer,DELIM)) {
        *token = '\0';
    }
    *h2 = atof(buffer);
    printf("%s,h2=%f\n",buffer,*h2);

    if (fgets(buffer,MAX_REA_NAME,ifp) == NULL) {
        error_msg_fatal("missing lpc info?\n");
    }
    if (token=strpbrk(buffer,DELIM)) {
        *token = '\0';
    }
    *lpc = atof(buffer);
    printf("%s,lpc=%f\n",buffer,*lpc);

    if (fgets(buffer,MAX_REA_NAME,ifp) == NULL) {
        error_msg_fatal("missing gpc info?\n");
    }
    if (token=strpbrk(buffer,DELIM)) {
        *token = '\0';
    }
    *gpc = atof(buffer);
    printf("%s,gpc=%f\n",buffer,*gpc);

    fclose(ifp);

#ifdef DEBUG
    error_msg_warning("hmt_fix_parameters() :: end\n");
#endif
}

