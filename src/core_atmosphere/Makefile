.SUFFIXES: .F .o

#PHYSICS=-DDO_PHYSICS
PHYSICS=-DDO_PHYSICS_GFS
#PHYSICS=

OBJS = mpas_atm_mpas_core.o \
       mpas_atm_interp_diagnostics.o

all: physcore dycore atmcore

core_reg:
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) Registry.xml > Registry_processed.xml

ifeq ($(PHYSICS),-DDO_PHYSICS)
physcore:
	( cd physics; $(MAKE) all )
	( mkdir libphys; cd libphys; ar -x ../physics/libphys.a )
	( cd ../..; ln -sf ./src/core_atmosphere/physics/physics_wrf/files/*TBL .)
	( cd ../..; ln -sf ./src/core_atmosphere/physics/physics_wrf/files/*DATA* .)
	( cd ../..; ln -sf ./src/core_atmosphere/physics/physics_wrf/files/*DBL .)
else
ifeq ($(PHYSICS),-DDO_PHYSICS_GFS)
physcore:
	( cd physics_gfs; $(MAKE) all )
	( mkdir libphys; cd libphys; ar -x ../physics_gfs/libphys.a )
endif
endif

dycore:
	( cd dynamics; $(MAKE) all PHYSICS="$(PHYSICS)" )

atmcore: physcore dycore $(OBJS)
	ar -ru libdycore.a $(OBJS) dynamics/*.o libphys/*.o

mpas_atm_mpas_core.o: dycore mpas_atm_interp_diagnostics.o

ifeq ($(PHYSICS),-DDO_PHYSICS)
clean:
	( cd physics; $(MAKE) clean )
	( cd dynamics; $(MAKE) clean )
	( cd ../..; rm -f *TBL )
	( cd ../..; rm -f *DATA* )
	$(RM) -r libphys
	$(RM) *.o *.mod *.f90 libdycore.a
	$(RM) Registry_processed.xml
	@# Certain systems with intel compilers generate *.i files
	@# This removes them during the clean process
	$(RM) *.i

.F.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(PHYSICS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../framework -I../operators -I./physics -I./dynamics -I./physics/physics_wrf -I../external/esmf_time_f90
else
	$(FC) $(CPPFLAGS) $(PHYSICS) $(FFLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES) -I../framework -I../operators -I./physics -I./dynamics -I./physics/physics_wrf -I../external/esmf_time_f90
endif

ifeq ($(PHYSICS),-DDO_PHYSICS_GFS)
clean:
	( cd physics_gfs; $(MAKE) clean )
	( cd dynamics; $(MAKE) clean )
	$(RM) -r libphys
	$(RM) *.o *.mod *.f90 libdycore.a

.F.o:
	$(RM) $@ $*.mod
	$(CPP) $(CPPFLAGS) $(PHYSICS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../framework -I../operators -I./physics_gfs -I./dynamics -I./physics_gfs/physics_gfs_lib -I../external/esmf_time_f90
endif

endif

