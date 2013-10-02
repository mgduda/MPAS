.SUFFIXES: .F .o

ifeq ($(CORE),atmosphere)
COREDEF = -Dmpas
endif

dummy:
	echo "****** make physics_gfs ******"

all: core_physics_gfs_lib core_physics_gfs_src

core_physics_gfs_lib:
	(cd physics_gfs_lib; make all COREDEF="$(COREDEF)")

core_physics_gfs_src:
	(cd physics_gfs_src; make all COREDEF="$(COREDEF)")

clean:
	$(RM) *.o *.mod *.f90 libphys.a
	( cd physics_gfs_lib; make clean )
	( cd physics_gfs_src; make clean )

.F.o:
	$(RM) $@ $*.mod
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../../framework -I../../operators -I./physics_gfs_lib -I./physics_gfs_src -I../../external/esmf_time_f90


