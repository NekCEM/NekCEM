?FC?
?FFLAGS?
?CC?
?CFLAGS?
?LD?
?LDFLAGS?
?NEKBASE?
?JLBASE?
?APP?
?USR?
?USR_LFLAGS?
?USR_IFLAGS?

BINNAME=nekcem
APPSRC = $(shell find $(NEKBASE) -name cem_$(APP)*.F)
APPOBJ = $(patsubst $(NEKBASE)/%.F, obj/%.o, $(APPSRC))
FSRC = $(shell find $(NEKBASE) -name *.F | grep -v 'nek5_comm' \
					 | grep -v 'maxwell' \
					 | grep -v 'drift' \
					 | grep -v 'schrod')
FOBJ = $(patsubst $(NEKBASE)/%.F, obj/%.o, $(FSRC))
CSRC = $(shell find $(NEKBASE) -maxdepth 1 -name *.c)
COBJ = $(patsubst $(NEKBASE)/%.c, obj/%.o, $(CSRC))
COMMSRC = $(NEKBASE)/nek5_comm_mpi.F
COMMOBJ = $(patsubst $(NEKBASE)/%.F, obj/%.o, $(COMMSRC))
# Unfortunately, for now there doesn't seem to be a sensible way to
# specify which files in jl we want to compile short of listing them
# all.
JLFILES = gs.c sarray_transfer.c sarray_sort.c gs_local.c crystal.c	\
comm.c tensor.c fail.c fcrystal.c findpts.c findpts_local.c obbox.c	\
poly.c lob_bnd.c findpts_el_3.c findpts_el_2.c sparse_cholesky.c	\
xxt.c fcrs.c sort.c
JLSRC = $(patsubst %, $(JLBASE)/%, $(JLFILES))
JLOBJ = $(patsubst $(JLBASE)/%.c, obj/%.o, $(JLSRC))
CASENAME = $(basename $(USR))
USROBJ = obj/${CASENAME}.o

obj/%.o: $(NEKBASE)/%.F SIZE
	$(FC) -c $(FFLAGS) $< -o $@

obj/%.o: $(NEKBASE)/%.c
	$(CC) -c $(CFLAGS) $(USR_IFLAGS) $< -o $@

obj/%.o: $(JLBASE)/%.c
	$(CC) -c $(CFLAGS) $< -o $@

all: usr jl $(APPOBJ) $(FOBJ) $(COBJ) $(COMMOBJ)
	$(LD) -o nekcem $(APPOBJ) $(FOBJ) $(COBJ) $(JLOBJ) $(USROBJ) \
			$(COMMOBJ) $(LDFLAGS) $(USR_LFLAGS)
	@if test -f ${BINNAME}; then \
        echo ""; \
        echo "$(NEK_WARN)"; \
        echo ""; \
	echo "#############################################################"; \
	echo "#                  Compilation successful!                  #"; \
	echo "#############################################################"; \
	size ${BINNAME}; \
	echo ""; \
	else \
	echo -e "\033[1;31;38m" "ERROR: Compilation failed!"; \
	echo -e "\033[0m"; \
	fi

usr: ${USROBJ}

${USROBJ}: ./SIZE ${USR}
	@env CASENAME=$(CASENAME) $(NEKBASE)/mkuserfile
	$(FC) -c ${CASENAME}.F -o $(USROBJ) $(FFLAGS) -I$(NEKBASE)

jl: dir $(JLOBJ)

dir:
	@mkdir -p obj vtk

clean:
	@rm -rf compiler.out nekcem obj box.tmp
	@rm -f xxt_map.rea SESSION.NAME
	@rm -f *.cobaltlog *.error logfile
	@find . -type f -name "*~" -exec rm {} \;

cleanall:
	@rm -rf compiler.out nekcem vtk obj box.tmp
	@rm -f *.output xxt_map.rea SESSION.NAME
	@rm -f *.cobaltlog *.error logfile
	@find . -type f -name "*~" -exec rm {} \;
