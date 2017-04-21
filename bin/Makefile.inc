?FC?
?FFLAGS?
?CC?
?CFLAGS?
?LD?
?LDFLAGS?
?NEKBASE?
?JLBASE?

FSRC = $(shell find $(NEKBASE) -name *.F | grep -v 'nek5_comm' \
					 | grep -v 'mxm')
FOBJ = $(patsubst $(NEKBASE)/%.F, obj/%.o, $(FSRC))
CSRC = $(shell find $(NEKBASE) -maxdepth 1 -name *.c)
COBJ = $(patsubst $(NEKBASE)/%.c, obj/%.o, $(CSRC))
COMMSRC = $(NEKBASE)/nek5_comm_mpi.F
COMMOBJ = $(patsubst $(NEKBASE)/%.F, obj/%.o, $(COMMSRC))
MXMSRC = $(NEKBASE)/nek5_mxm_std.F $(NEKBASE)/nek5_mxm_wrapper.F
MXMOBJ = $(patsubst $(NEKBASE)/%.F, obj/%.o, $(MXMSRC))
# Unfortunately, for now there doesn't seem to be a sensible way to
# specify which files in jl we want to compile short of listing them
# all.
JLFILES = gs.c sarray_transfer.c sarray_sort.c gs_local.c crystal.c	\
comm.c tensor.c fail.c fcrystal.c findpts.c findpts_local.c obbox.c	\
poly.c lob_bnd.c findpts_el_3.c findpts_el_2.c sparse_cholesky.c	\
xxt.c fcrs.c sort.c
JLSRC = $(patsubst %, $(JLBASE)/%, $(JLFILES))
JLOBJ = $(patsubst $(JLBASE)/%.c, obj/%.o, $(JLSRC))
USR = $(wildcard *.usr)
USROBJ = obj/subuser.o


obj/%.o : $(NEKBASE)/%.F
	$(FC) -c $(FFLAGS) $< -o $@

obj/%.o : $(NEKBASE)/%.c
	$(CC) -c $(CFLAGS) $< -o $@

obj/%.o : $(JLBASE)/%.c
	$(CC) -c $(CFLAGS) $< -o $@

all : usr jl $(FOBJ) $(COBJ) $(COMMOBJ) $(MXMOBJ)
	$(LD) -o nekcem $(FOBJ) $(COBJ) $(JLOBJ) $(USROBJ) $(COMMOBJ) \
	                $(MXMOBJ) $(LDFLAGS)

usr : dir $(USR)
	cp $(USR) obj/subuser.F
	$(FC) -c obj/subuser.F -o $(USROBJ) $(FFLAGS) -I$(NEKBASE)
	rm obj/subuser.F

jl : CFLAGS += -DAMG_DUMP -DNEW_GS_LOOPS
jl : dir $(JLOBJ)

dir :
	@mkdir -p obj vtk

clean :
	rm -rf compiler.out nekcem vtk wake obj
	rm -f *.output xxt_map.rea SESSION.NAME
	rm -f *.cobaltlog *.error logfile
	find . -type f -name "*~" -exec rm {} \;
