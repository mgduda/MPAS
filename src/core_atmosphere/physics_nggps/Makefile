.SUFFIXES: .F .F90 .o

ifeq ($(CORE),atmosphere)
COREDEF = -Dmpas
endif

dummy:
	echo "****** make physics_nggps ******"

OBJS = \
	mpas_atmphys_nggps_cloud_properties.o \
	mpas_atmphys_nggps_constants.o \
	mpas_atmphys_nggps_control.o \
	mpas_atmphys_nggps_date_time.o \
	mpas_atmphys_nggps_diagnostics.o \
	mpas_atmphys_nggps_driver.o \
	mpas_atmphys_nggps_driver_out.o \
	mpas_atmphys_nggps_dynamic_parameters.o \
	mpas_atmphys_nggps_init.o \
	mpas_atmphys_nggps_interface_fields.o \
	mpas_atmphys_nggps_manager.o \
	mpas_atmphys_nggps_model_parameters.o \
	mpas_atmphys_nggps_o3climatology.o \
	mpas_atmphys_nggps_radiation_tendencies.o \
	mpas_atmphys_nggps_sfc_properties.o \
	mpas_atmphys_nggps_state_fields_in.o \
	mpas_atmphys_nggps_state_fields_out.o \
	mpas_atmphys_nggps_tbd_ddt.o \
	mpas_atmphys_nggps_todynamics.o \
	mpas_atmphys_nggps_update.o \
        mpas_atmphys_nggps_update_surface.o \
	mpas_atmphys_nggps_vars.o \
	mpas_atmphys_date_time.o \
	mpas_atmphys_utilities.o \
	nuopc_physics.o

all: core_physics_nggps_lib core_physics_nggps_src core_physics_nggps

core_physics_nggps_lib:
	(cd physics_nggps_lib; make all COREDEF="$(COREDEF)")

core_physics_nggps_src:
	(cd physics_nggps_src; make all COREDEF="$(COREDEF)")

core_physics_nggps: $(OBJS)
	ar -ru libphys.a $(OBJS)

# DEPENDENCIES:
mpas_atmphys_nggps_cloud_properties.o: \
	mpas_atmphys_nggps_vars.o \
	nuopc_physics.o

mpas_atmphys_nggps_date_time.o: \
	mpas_atmphys_date_time.o

mpas_atmphys_nggps_diagnostics.o: \
	mpas_atmphys_nggps_vars.o \
	nuopc_physics.o \
	./physics_nggps_src/machine.o \
	./physics_nggps_src/radlw_param.o \
	./physics_nggps_src/radsw_param.o

mpas_atmphys_nggps_driver.o: \
	mpas_atmphys_nggps_cloud_properties.o \
	mpas_atmphys_nggps_diagnostics.o \
	mpas_atmphys_nggps_driver_out.o \
	mpas_atmphys_nggps_dynamic_parameters.o \
	mpas_atmphys_nggps_interface_fields.o \
        mpas_atmphys_nggps_radiation_tendencies.o \
	mpas_atmphys_nggps_sfc_properties.o \
	mpas_atmphys_nggps_state_fields_in.o \
	mpas_atmphys_nggps_state_fields_out.o \
	mpas_atmphys_nggps_tbd_ddt.o \
	mpas_atmphys_nggps_vars.o \
	nuopc_physics.o

mpas_atmphys_nggps_driver_out.o: \
	mpas_atmphys_nggps_vars.o \
	nuopc_physics.o \
	./physics_nggps_src/machine.o \
	./physics_nggps_src/radlw_param.o \
	./physics_nggps_src/radsw_param.o

mpas_atmphys_nggps_driver_radiation.o: \
	mpas_atmphys_nggps_dynamic_parameters.o \
	mpas_atmphys_nggps_vars.o \
	nuopc_physics.o

mpas_atmphys_nggps_dynamic_parameters.o: \
	mpas_atmphys_nggps_manager.o \
	mpas_atmphys_nggps_vars.o \
	nuopc_physics.o \
	./physics_nggps_src/grrad.o \
	./physics_nggps_src/mersenne_twister.o

mpas_atmphys_nggps_init.o: \
	mpas_atmphys_nggps_model_parameters.o \
	mpas_atmphys_nggps_o3climatology.o \
	mpas_atmphys_nggps_vars.o \
	nuopc_physics.o  \
	./physics_nggps_src/machine.o

mpas_atmphys_interface_fields.o: \
	nuopc_physics.o \
	./physics_nggps_src/machine.o \
	./physics_nggps_src/radlw_param.o \
	./physics_nggps_src/radsw_param.o

mpas_atmphys_nggps_manager.o : \
	mpas_atmphys_nggps_o3climatology.o \
	mpas_atmphys_nggps_update.o \
	mpas_atmphys_nggps_update_surface.o \
	mpas_atmphys_nggps_vars.o \
	mpas_atmphys_utilities.o

mpas_nggps_model_parameters.o: \
	./physics_nggps_src/machine.o \
	./physics_nggps_src/ozne_def.o

mpas_atmphys_nggps_o3climatology.o: \
	mpas_atmphys_nggps_constants.o \
	mpas_atmphys_nggps_vars.o

mpas_atmphys_nggps_radiation_tendencies.o: \
	nuopc_physics.o

mpas_atmphys_nggps_sfc_properties.o: \
	mpas_atmphys_nggps_vars.o \
	nuopc_physics.o \
	./physics_nggps_src/machine.o

mpas_atmphys_nggps_state_fields_in.o: \
	mpas_atmphys_nggps_constants.o \
	mpas_atmphys_nggps_vars.o \
	nuopc_physics.o \
	./physics_nggps_src/machine.o

mpas_atmphys_nggps_state_fields_out.o: \
	mpas_atmphys_nggps_vars.o \
	nuopc_physics.o

mpas_atmphys_nggps_tbd_ddt.o: \
	mpas_atmphys_nggps_vars.o \
	nuopc_physics.o \
	./physics_nggps_src/machine.o

mpas_atmphys_nggps_todynamics.o: \
	mpas_atmphys_nggps_constants.o

mpas_atmphys_nggps_update_surface.o: \
	mpas_atmphys_nggps_date_time.o \
	mpas_atmphys_date_time.o

mpas_atmphys_nggps_vars.o: \
	nuopc_physics.o \
	./physics_nggps_src/machine.o \
	./physics_nggps_src/funcphys.o \
	./physics_nggps_src/mersenne_twister.o \
	./physics_nggps_src/physcons.o \
	./physics_nggps_src/physparam.o \
	./physics_nggps_src/ozne_def.o \
	./physics_nggps_src/rascnvv2.o

nuopc_physics.o: \
	./physics_nggps_src/grrad.o \
	./physics_nggps_src/machine.o \
	./physics_nggps_src/physcons.o \
	./physics_nggps_src/radlw_param.o \
	./physics_nggps_src/radsw_param.o

clean:
	$(RM) *.o *.mod *.f90 libphys.a
	( cd physics_nggps_lib; make clean )
	( cd physics_nggps_src; make clean )

.F.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../../framework -I./physics_nggps_lib -I./physics_nggps_src -I../../external/esmf_time_f90
else
	$(FC) $(CPPFLAGS) $(COREDEF) $(FFLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES) -I../../framework -I./physics_nggps_lib -I./physics_nggps_src -I../../external/esmf_time_f90
endif

.F90.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../../framework -I./physics_nggps_lib -I./physics_nggps_src -I../../external/esmf_time_f90
else
	$(FC) $(CPPFLAGS) $(COREDEF) $(FFLAGS) -c $*.F90 $(CPPINCLUDES) $(FCINCLUDES) -I../../framework -I./physics_nggps_lib -I./physics_nggps_src -I../../external/esmf_time_f90
endif

