###############################################################################
L	= 
S       = $(BASEDIR)/src
NETLIB	= $(BASEDIR)/libs/netlib
D	= $(BASEDIR)/libs/libtfs
R	= $(BASEDIR)/libs/libtfs/rsb
################################################################################
#
# core nekton files ...
#
#
CORE = coef.o connect1.o connect2.o drive.o \
bdry.o byte.o vtkbin.o dssum.o edgec.o genxyz.o \
prepost.o speclib.o subs1.o subs2.o subuser.o \
convert.o map2.o mat1.o mat2std.o \
genbox.o ic.o movwin.o\
grad.o cem_dg.o cem_opt.o bessel.o randgen.o pml.o dumpvtk.o \
main.o dnchbv.o outmesh.o poisson2d.o beam.o wake.o courant.o\
#
# pff overlapping schwaz 
#
OS =
#outmesh.o 
#fastsa.o 
#
# gather-scatter package 
# 
GS = bit_mask.o bss_malloc.o ivec.o error.o comm.o gs.o queue.o hmt_perf.o flush_hack.o 
#
# coarse grid solver
XXT = xxt.o xxt_tools.o tools.o
RSB = pythag.o tqli.o smbwr.o  sparse_matrix.o stack.o adj_list.o
#
#
#
NOBJS = $(CORE) $(OS) $(GS) $(XXT) $(RSB) mxm44f2.o $(COMM_BACKEND).o \
	$(LAPACK_EXTRA_OBJS)

L0 =  $(OPT_OPTIONS) -O0
L2 =  $(OPT_OPTIONS) -O2 
L3 =  $(OPT_OPTIONS) -O3
L4 =  $(OPT_OPTIONS) -O4
HMT   = #-Dr8 -DDEBUG -DMPISRC -DSINGLE_DIR # -DBLAS

###
FL0   = $(L0) $(HMT) $(P)
FL2   = $(L2) $(HMT) $(P)
FL3   = $(L3) $(HMT) $(P)
FL4   = $(L4) $(HMT) $(P)

###
DRL0    = $(FL0) # -Mconcur -Mreentrant
DFL0    = $(FL0) # -Mconcur 
DFL2    = $(FL2) # -Mconcur
DFL3    = $(FL3) # -Mconcur 
DFL4    = $(FL4) # -Mconcur -Minfo=loop -Mneginfo=loop

RFL0    = $(FL0) #-Mreentrant
RFL2    = $(FL2) #-Mreentrant
RFL3    = $(FL3) #-Mreentrant
RFL4    = $(FL4) #-Mreentrant

DRFL0   = $(FL0) # -Mconcur -Mreentrant
DRFL2   = $(FL3) # -Mconcur -Mreentrant
DRFL3   = $(FL3) # -Mconcur -Mreentrant
DRFL4   = $(FL4) # -Mconcur -Mreentrant

cFL0   = $(L0) $(HMT)
cFL2   = $(L2) $(HMT)
cFL3   = $(L3) $(HMT)
cFL4   = $(L4) $(HMT)

###

lFLAGS = $(LAPACK) $(MPI_LIBS) $(EXTRALIBS)

################################################################################

all : nekcem  print

nekcem:	 $(NOBJS) 
	 $(LD) $(OPT_OPTIONS) -o nekcem $G $(NOBJS) $(lFLAGS) 
# -Bstatic

print: nekcem  
	@if which size &> /dev/null ; then size nekcem ; fi
	@echo "I am done" 

clean:
	rm *.o nekcem  


################################################################################
# fortran source files
################################################################################
drive.o		:$S/drive.F;		$(F77) -c $(FL2) $<
main.o		:$S/main.c;		$(CC) -c $<
prepost.o	:$S/prepost.F;		$(F77) -c $(FL2) $<
connect1.o	:$S/connect1.F;		$(F77) -c $(FL2) $<
connect2.o	:$S/connect2.F;		$(F77) -c $(FL2) $<
edgec.o		:$S/edgec.F;		$(F77) -c $(FL2) $<
ginterface.o	:$S/ginterface.F;	$(F77) -c $(FL2) $<
genxyz.o	:$S/genxyz.F;		$(F77) -c $(FL2) $<
genbox.o	:$S/genbox.F;		$(F77) -c $(FL2) $<
dumpvtk.o	:$S/dumpvtk.F;		$(F77) -c $(FL2) $<
poisson2d.o	:$S/poisson2d.F;	$(F77) -c $(FL2) $<
beam.o	        :$S/beam.F;	        $(F77) -c $(FL2) $<
wake.o	        :$S/wake.F;	        $(F77) -c $(FL2) $<
movwin.o        :$S/movwin.F;	        $(F77) -c $(FL2) $<
courant.o       :$S/courant.F;	        $(F77) -c $(FL2) $<
subs1.o		:$S/subs1.F;		$(F77) -c $(FL2) $<
subs2.o		:$S/subs2.F;		$(F77) -c $(FL2) $<
turb.o		:$S/turb.F;		$(F77) -c $(FL2) $<
map2.o		:$S/map2.F;		$(F77) -c $(FL2) $<
mvmesh.o	:$S/mvmesh.F;		$(F77) -c $(FL2) $<
reverf.o	:$S/reverf.F;		$(F77) -c $(FL2) $<
################################################################################
bdry.o		:$S/bdry.F;		$(F77) -c $(FL2) $<
commcube.o	:$S/commcube.F;		$(F77) -c $(FL2) $<
comm_seq.o	:$S/comm_seq.F;		$(F77) -c $(FL2) $<
comm_mpi.o	:$S/comm_mpi.F;		$(F77) -c $(FL2) $<
subuser.o	:$S/subuser.F;		$(F77) -c $(FL2) $<
ic.o		:$S/ic.F;		$(F77) -c $(FL2) $<
conduct.o	:$S/conduct.F;		$(F77) -c $(FL2) $<
col2a.o		:$S/col2a.s;		$(F77) -c $(FL2) $<
col3a.o		:$S/col3a.s;		$(F77) -c $(FL2) $<
add2a.o		:$S/add2a.s;		$(F77) -c $(FL2) $<
add3a.o		:$S/add3a.s;		$(F77) -c $(FL2) $<
convert.o	:$S/convert.F;		$(F77) -c $(FL2) $<
outmesh.o	:$S/outmesh.F;		$(F77) -c $(FL2) $<
navier0.o	:$S/navier0.F;		$(F77) -c $(FL2) $<
navier2.o	:$S/navier2.F;		$(F77) -c $(FL2) $<
navier3.o	:$S/navier3.F;		$(F77) -c $(FL2) $<
navier4.o	:$S/navier4.F;		$(F77) -c $(FL2) $<
navier5.o	:$S/navier5.F;		$(F77) -c $(FL2) $<
navier6.o	:$S/navier6.F;		$(F77) -c $(FL2) $<
navier7.o	:$S/navier7.F;		$(F77) -c $(FL2) $<
navier8.o	:$S/navier8.F;		$(F77) -c $(FL2) $<
speclib.o	:$S/speclib.F;		$(F77) -c $(FL2) $<
calcz.o		:$S/calcz.F;		$(F77) -c $(FL2) $<
fast3d.o	:$S/fast3d.F;		$(F77) -c $(FL2) $<
util_pc.o       :$S/util_pc.F;		$(F77) -c $(FL2) $<
qqt.o           :$S/qqt.F;		$(F77) -c $(FL2) $<
hmholtz_pc.o    :$S/hmholtz_pc.F;	$(F77) -c $(FL2) $<
operator.o      :$S/operator.F;		$(F77) -c $(FL2) $<
setup_pc1.o     :$S/setup_pc1.F;	$(F77) -c $(FL2) $<
setup_pc2.o     :$S/setup_pc2.F;	$(F77) -c $(FL2) $<
################################################################################
dssum.o		:$S/dssum.F;		$(F77) -c $(FL3) $<
dssmr.o		:$S/dssmr.F;		$(F77) -c $(FL3) $<
dssm0.o		:$S/dssm0.F;		$(F77) -c $(FL3) $<
eigsolv.o	:$S/eigsolv.F;		$(F77) -c $(FL3) $<
hmh_dg.o	:$S/hmh_dg.F;		$(F77) -c $(FL3) $<
cem_dg.o	:$S/cem_dg.F;		$(F77) -c $(FL3) $<
cem_opt.o	:$S/cem_opt.F;		$(F77) -c $(FL3) $<
randgen.o	:$S/randgen.F;		$(F77) -c $(FL3) $<
pml.o    	:$S/pml.F;		$(F77) -c $(FL3) $<
bessel.o	:$S/bessel.F;		$(F77) -c $(FL3) $<
dnchbv.o	:$S/dnchbv.F;		$(F77) -c $(FL3) $<
grad.o		:$S/grad.F;		$(F77) -c $(FL3) $<
gfdm_par.o	:$S/gfdm_par.F;		$(F77) -c $(FL3) $<
gfdm_solve.o	:$S/gfdm_solve.F;	$(F77) -c $(FL3) $<
gfdm_op.o	:$S/gfdm_op.F;		$(F77) -c $(FL3) $<
gauss.o		:$S/gauss.F;		$(F77) -c $(FL3) $<
planx.o		:$S/planx.F;		$(F77) -c $(FL3) $<
ssolv.o		:$S/ssolv.F;		$(F77) -c $(FL3) $<
mperre.o	:$S/mperre.F;		$(F77) -c $(FL3) $<
per2.o          :$S/per2.F;             $(F77) -c $(FL3) $<
################################################################################
coef.o		:$S/coef.F;		$(F77) -c $(FL4) $<
hsmg.o		:$S/hsmg.F;		$(F77) -c $(FL4) $<
uzawa_gmres.o	:$S/uzawa_gmres.F;	$(F77) -c $(FL4) $<
mat1.o		:$S/mat1.F;		$(F77) -c $(FL4) $<
mat2.o		:$S/mat2.F;		$(F77) -c $(FL4) $<
mat2std.o	:$S/mat2std.F;		$(F77) -c $(FL4) $<
xxt_tools.o     :$S/xxt_tools.F;        $(F77) -c $(FL4) $<

################################################################################
# NETLIB
################################################################################
blas.o          :$(NETLIB)/blas.f;      $(F77) -c $(FL4) $<
lapack.o        :$(NETLIB)/lapack.f;    $(F77) -c $(FL4) $<
arpack.o        :$(NETLIB)/arpack.f;    $(F77) -c $(FL4) $<


################################################################################
#reentrant modules
################################################################################
mxm.o		:$S/mxm.F;		$(F77) -c $(RFL4) $<
mxm44.o		:$S/mxm44.F;		$(F77) -c $(RFL4) $<
mxm44f2.o	:$S/mxm44f2.F;		$(F77) -c $(RFL4) $<
hmholtza.o	:$S/hmholtza.F;		$(F77) -c $(RFL4) $<
navier1a.o	:$S/navier1a.F;		$(F77) -c $(RFL4) $<
fastsa.o	:$S/fastsa.F;		$(F77) -c $(RFL4) $<
################################################################################


################################################################################
#dual mode modules
################################################################################
hmholtz.o	:$S/hmholtz.F;		$(F77) -c $(DFL4) $<
navier1.o	:$S/navier1.F;		$(F77) -c $(DFL4) $<
fasts.o		:$S/fasts.F;		$(F77) -c $(DFL4) $<
fasts2.o	:$S/fasts2.F;		$(F77) -c $(RFL4) $<
#subuser.o	:$S/subuser.F;		$(F77) -c $(FL2) $<
################################################################################
################################################################################
# c source files
################################################################################
################################################################################
revert.o        :$S/revert.c;           $(CC)  -c $(cFL3) $<
byte.o          :$S/byte.c;             $(CC)  -c $(cFL3) $<
vtkbin.o    :$S/vtkbin.c;       $(CC)  -c $(cFL3) $<
hmt_perf.o      :$S/hmt_perf.c;         $(CC)  -c -I$D $(cFL3) -I$L/include $<
bit_mask.o      :$D/bit_mask.c;         $(CC)  -c $(cFL3) $<
bss_malloc.o    :$D/bss_malloc.c;       $(CC)  -c $(cFL3) -I$L/include $<
error.o         :$D/error.c;            $(CC)  -c $(cFL3) -I$L/include $<
queue.o         :$D/queue.c;            $(CC)  -c $(cFL3) $<
stack.o         :$D/stack.c;            $(CC)  -c $(cFL3) $<
sparse_matrix.o :$R/sparse_matrix.c;	$(CC)  -c -I$D $(cFL3) -I$L/include $<
tqli.o          :$R/tqli.c;         	$(CC)  -c -I$D $(cFL3) -I$L/include $<
pythag.o        :$R/pythag.c;       	$(CC)  -c -I$D $(cFL3) -I$L/include $<
adj_list.o      :$R/adj_list.c;     	$(CC)  -c -I$D $(cFL3) -I$L/include $<
smbwr.o         :$R/smbwr.c;        	$(CC)  -c -I$D $(cFL3) -I$L/include $<
################################################################################
xxt.o           :$D/xxt.c;              $(CC)  -c $(cFL4) -I$L/include $<
tools.o         :$D/tools.c;            $(CC)  -c $(cFL4) -I$L/include $<
ivec.o          :$D/ivec.c;             $(CC)  -c $(cFL4) -I$L/include $<
gs.o            :$D/gs.c;               $(CC)  -c $(cFL4) -I$L/include $<
comm.o          :$D/comm.c;             $(CC)  -c $(cFL4) -I$L/include $<
################################################################################
flush_hack.o    :$S/flush_hack.c;       $(CC)  -c -I$D -I$L/include $(cFL4) $<
mpi_hack.o      :$S/mpi_hack.c;         $(CC)  -c -I$D -I$L/include $(cFL4) $<
