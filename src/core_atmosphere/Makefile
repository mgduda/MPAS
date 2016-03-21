.SUFFIXES: .F .o

#PHYSICS=-DDO_PHYSICS
PHYSICS=-DDO_PHYSICS_NGGPS
#PHYSICS=

OBJS = mpas_atm_core.o \
       mpas_atm_core_interface.o \
       mpas_atm_interp_diagnostics.o \
       mpas_atm_dimensions.o

all: physcore dycore atmcore

core_reg:
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) Registry.xml > Registry_processed.xml

core_input_gen:
	if [ ! -e default_inputs ]; then  mkdir default_inputs; fi
	( cd default_inputs; $(NL_GEN) ../Registry_processed.xml namelist.atmosphere in_defaults=true )
	( cd default_inputs; $(ST_GEN) ../Registry_processed.xml streams.atmosphere stream_list.atmosphere. listed )

gen_includes:
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) Registry.xml > Registry_processed.xml
	(if [ ! -d inc ]; then mkdir -p inc; fi) # To generate *.inc files
	(cd inc; $(REG_PARSE) < ../Registry_processed.xml )

post_build:
	if [ ! -e $(ROOT_DIR)/default_inputs ]; then mkdir $(ROOT_DIR)/default_inputs; fi
	cp default_inputs/* $(ROOT_DIR)/default_inputs/.
	( cd $(ROOT_DIR)/default_inputs; for FILE in `ls -1`; do if [ ! -e ../$$FILE ]; then cp $$FILE ../.; fi; done )

ifeq ($(PHYSICS),-DDO_PHYSICS)
physcore:
	( cd physics; $(MAKE) all )
	( mkdir libphys; cd libphys; ar -x ../physics/libphys.a )
	( cd ../..; ln -sf ./src/core_atmosphere/physics/physics_wrf/files/*TBL .)
	( cd ../..; ln -sf ./src/core_atmosphere/physics/physics_wrf/files/*DATA* .)
else
ifeq ($(PHYSICS),-DDO_PHYSICS_NGGPS)
physcore: mpas_atm_dimensions.o
	( cd physics_nggps; $(MAKE) all )
	( mkdir libphys; cd libphys; ar -x ../physics_nggps/libphys.a )
	( cd ../..; ln -sf ./src/core_atmosphere/physics_nggps/physics_nggps_files/*.f77 .)
	( cd ../..; ln -sf ./src/core_atmosphere/physics_nggps/physics_nggps_files/global_climaeropac_global.txt aerosol.dat)
	( cd ../..; ln -sf ./src/core_atmosphere/physics_nggps/physics_nggps_files/global_sfc_emissivity_idx.txt sfc_emissivity_idx.txt)
	( cd ../..; ln -sf ./src/core_atmosphere/physics_nggps/physics_nggps_files/global_solarconstant_noaa_an.txt solarconstant_noaa_an.txt)
	( cd ../..; ln -sf ./src/core_atmosphere/physics_nggps/physics_nggps_files/fix_co2_proj/global_co2historicaldata_2013.txt co2historicaldata_2013.txt)
	( cd ../..; ln -sf ./src/core_atmosphere/physics_nggps/physics_nggps_files/fix_co2_proj/global_co2historicaldata_2014.txt co2historicaldata_2014.txt)
	( cd ../..; ln -sf ./src/core_atmosphere/physics_nggps/physics_nggps_files/fix_co2_proj/global_co2historicaldata_2015.txt co2historicaldata_2015.txt)
	( cd ../..; ln -sf ./src/core_atmosphere/physics_nggps/physics_nggps_files/fix_co2_proj/global_co2historicaldata_2016.txt co2historicaldata_2016.txt)
endif
endif
	( cd ../..; ln -sf ./src/core_atmosphere/avec/avec_timers_config.txt avec_timers_config.txt )

dycore: mpas_atm_dimensions.o
	( cd dynamics; $(MAKE) all PHYSICS="$(PHYSICS)" )

atmcore: physcore dycore $(OBJS)
	ar -ru libdycore.a $(OBJS) dynamics/*.o libphys/*.o

mpas_atm_core_interface.o: mpas_atm_core.o

mpas_atm_core.o: dycore mpas_atm_interp_diagnostics.o

mpas_atm_dimensions.o:

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
	$(RM) -r default_inputs

.F.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(PHYSICS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../framework -I../operators -I./physics -I./dynamics -I./physics/physics_wrf -I../external/esmf_time_f90
else
	$(FC) $(CPPFLAGS) $(PHYSICS) $(FFLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES) -I../framework -I../operators -I./physics -I./dynamics -I./physics/physics_wrf -I../external/esmf_time_f90
endif

else

ifeq ($(PHYSICS),-DDO_PHYSICS_NGGPS)
clean:
	( cd physics_nggps; $(MAKE) clean )
	( cd dynamics; $(MAKE) clean )
	( cd ../..; rm -f *.f77 )
	( cd ../..; rm -f *.dat )
	( cd ../..; rm -f *.txt )
	$(RM) -r libphys
	$(RM) *.o *.mod *.f90 libdycore.a

.F.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(PHYSICS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../framework -I../operators -I./physics_nggps -I./dynamics -I./physics_nggps/physics_nggps_src -I../external/esmf_time_f90
else
	$(FC) $(CPPFLAGS) $(PHYSICS) $(FFLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES) -I../framework -I../operators -I./physics_nggps -I./dynamics -I./physics_nggps/physics_nggps_src -I../external/esmf_time_f90
endif
endif
endif


