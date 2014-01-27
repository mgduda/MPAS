.SUFFIXES: .F .o

OBJS = mpas_atm_time_integration.o \
       mpas_atm_advection.o

all: $(OBJS)

mpas_atm_time_integration.o: 

mpas_atm_advection.o: 

clean:
	$(RM) *.o *.mod *.f90
	@# Certain systems with intel compilers generate *.i files
	@# This removes them during the clean process
	$(RM) *.i

.F.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
ifeq ($(PHYSICS),-DDO_PHYSICS)
	$(CPP) $(CPPFLAGS) $(PHYSICS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../../framework -I../../operators -I../physics -I../physics/physics_wrf -I../../external/esmf_time_f90
else
ifeq ($(PHYSICS),-DDO_PHYSICS_GFS)
	$(CPP) $(CPPFLAGS) $(PHYSICS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../../framework -I../../operators -I../physics_gfs -I../physics_gfs/physics_gfs_src -I../../external/esmf_time_f90
endif
endif

else
ifeq ($(PHYSICS),-DDO_PHYSICS)
	$(FC) $(CPPFLAGS) $(PHYSICS) $(FFLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES) -I../../framework -I../../operators -I../physics -I../physics/physics_wrf -I../../external/esmf_time_f90
else
ifeq ($(PHYSICS),-DDO_PHYSICS_GFS)
	$(FC) $(CPPFLAGS) $(PHYSICS) $(FFLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES) -I../../framework -I../../operators -I../physics_gfs -I../physics_gfs/physics_gfs_src -I../../external/esmf_time_f90
endif
endif
endif
