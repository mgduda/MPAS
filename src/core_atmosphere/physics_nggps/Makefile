.SUFFIXES: .F .F90 .o

ifeq ($(CORE),atmosphere)
COREDEF = -Dmpas
endif

dummy:
	echo "****** make physics_nggps ******"

OBJS = \
	nuopc_physics.o

all: core_physics_nggps_src core_physics_nggps

core_physics_nggps_src:
	(cd physics_nggps_src; make all COREDEF="$(COREDEF)")

core_physics_nggps: $(OBJS)
	ar -ru libphys.a $(OBJS)

# DEPENDENCIES:
nuopc_physics.o : \
	./physics_nggps_src/grrad.o \
	./physics_nggps_src/machine.o \
	./physics_nggps_src/physcons.o \
	./physics_nggps_src/radlw_param.o \
	./physics_nggps_src/radsw_param.o

clean:
	$(RM) *.o *.mod *.f90 libphys.a
	( cd physics_nggps_src; make clean )

.F.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../../framework -I./physics_nggps_src -I../../external/esmf_time_f90
else
	$(FC) $(CPPFLAGS) $(COREDEF) $(FFLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES) -I../../framework -I./physics_nggps_src -I../../external/esmf_time_f90
endif

.F90.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../../framework -I./physics_nggps_src -I../../external/esmf_time_f90
else
	$(FC) $(CPPFLAGS) $(COREDEF) $(FFLAGS) -c $*.F90 $(CPPINCLUDES) $(FCINCLUDES) -I../../framework -I./physics_nggps_src -I../../external/esmf_time_f90
endif

