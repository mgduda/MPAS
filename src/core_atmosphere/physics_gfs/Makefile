.SUFFIXES: .F .o

ifeq ($(CORE),atmosphere)
COREDEF = -Dmpas
endif

dummy:
	echo "****** make physics_gfs ******"

OBJS = \
	mpas_atmphys_date_time.o \
	mpas_atmphys_gfs_constants.o \
	mpas_atmphys_gfs_control.o \
	mpas_atmphys_gfs_date_time.o \
	mpas_atmphys_gfs_driver.o \
	mpas_atmphys_gfs_driver_gbphys.o \
	mpas_atmphys_gfs_driver_grrad.o \
	mpas_atmphys_gfs_grrad_init.o \
	mpas_atmphys_gfs_init.o \
	mpas_atmphys_gfs_interface.o \
	mpas_atmphys_gfs_manager.o \
	mpas_atmphys_gfs_o3climatology.o \
	mpas_atmphys_gfs_todynamics.o \
	mpas_atmphys_gfs_update.o \
	mpas_atmphys_gfs_update_surface.o \
	mpas_atmphys_gfs_vars.o \
	mpas_atmphys_utilities.o

all: core_physics_gfs_lib core_physics_gfs_src core_physics_gfs

core_physics_gfs: $(OBJS)
	ar -ru libphys.a $(OBJS)

core_physics_gfs_lib:
	(cd physics_gfs_lib; make all COREDEF="$(COREDEF)")

core_physics_gfs_src:
	(cd physics_gfs_src; make all COREDEF="$(COREDEF)")

#DEPENDENCIES:
mpas_atmphys_gfs_date_time.o: \
	mpas_atmphys_date_time.o

mpas_atmphys_gfs_driver.o: \
	mpas_atmphys_gfs_driver_gbphys.o \
	mpas_atmphys_gfs_driver_grrad.o \
	mpas_atmphys_gfs_interface.o \
	mpas_atmphys_gfs_manager.o \
	mpas_atmphys_gfs_vars.o

mpas_atmphys_gfs_driver_gbphys.o: \
	mpas_atmphys_gfs_interface.o \
	mpas_atmphys_gfs_manager.o \
	mpas_atmphys_gfs_vars.o \
	./physics_gfs_src/mersenne_twister.o \
	./physics_gfs_src/namelist_def.o

mpas_atmphys_gfs_driver_grrad.o: \
	mpas_atmphys_gfs_constants.o \
	mpas_atmphys_gfs_interface.o \
	mpas_atmphys_gfs_manager.o \
	mpas_atmphys_gfs_vars.o \
	mpas_atmphys_utilities.o \
	./physics_gfs_src/grrad.o \
	./physics_gfs_src/mersenne_twister.o \
	./physics_gfs_src/radiation_astronomy.o \
	./physics_gfs_src/radlw_param.o \
	./physics_gfs_src/radsw_param.o

mpas_atmphys_gfs_grrad_init.o: \
	mpas_atmphys_gfs_vars.o

mpas_atmphys_gfs_init.o: \
	mpas_atmphys_gfs_grrad_init.o \
	mpas_atmphys_gfs_o3climatology.o \
	mpas_atmphys_gfs_vars.o \
	./physics_gfs_src/funcphys.o \
	./physics_gfs_src/mersenne_twister.o \
	./physics_gfs_src/tracer_const_h-new.o

mpas_atmphys_gfs_interface.o: \
	mpas_atmphys_gfs_constants.o \
	mpas_atmphys_gfs_vars.o

mpas_atmphys_gfs_manager.o: \
	mpas_atmphys_gfs_o3climatology.o \
	mpas_atmphys_gfs_update_surface.o \
	mpas_atmphys_gfs_vars.o \
	mpas_atmphys_gfs_update.o \
	mpas_atmphys_utilities.o

mpas_atmphys_gfs_o3climatology.o: \
	mpas_atmphys_gfs_constants.o \
	mpas_atmphys_gfs_vars.o

mpas_atmphys_gfs_todynamics.o: \
	mpas_atmphys_gfs_constants.o

mpas_atmphys_gfs_update_surface.o: \
	mpas_atmphys_date_time.o \
	mpas_atmphys_gfs_date_time.o

mpas_atmphys_gfs_vars.o: \
        ./physics_gfs_src/date_def.o \
        ./physics_gfs_src/layout_grid_tracers.o \
        ./physics_gfs_src/mpi_def.o \
        ./physics_gfs_src/namelist_def.o \
        ./physics_gfs_src/namelist_soilveg.o \
	./physics_gfs_src/radlw_param.o \
	./physics_gfs_src/radsw_param.o \
        ./physics_gfs_src/resol_def.o

clean:
	$(RM) *.o *.mod *.f90 libphys.a
	( cd physics_gfs_lib; make clean )
	( cd physics_gfs_src; make clean )

.F.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../../framework -I../../operators -I./physics_gfs_lib -I./physics_gfs_src -I../../external/esmf_time_f90
else
	$(FC) $(CPPFLAGS) $(COREDEF) $(FFLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES) -I../../framework -I../../operators -I./physics_gfs_lib -I./physics_gfs_src -I../../external/esmf_time_f90
endif

