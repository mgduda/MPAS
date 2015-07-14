.SUFFIXES: .F .f .o

all: dummy physics_nggps_src

dummy:
	echo "****** compile physics_nggps_src ******"

OBJS = \
	calpreciptype.o \
	cnvc90.o \
	date_def.o \
	dcyc2.o \
	dcyc2.pre.rad.o \
	funcphys.o \
	gbphys.o \
	get_prs.o \
	gfs_phy_tracer_config.o \
	grrad.o \
	gscond.o \
	gscondp.o \
	gsmddrive.o \
	gwdc.o \
	gwdps.o \
	ideaca.o \
        iounitdef.o \
	lrgsclr.o \
	machine.o  \
	mersenne_twister.o \
	module_bfmicrophysics.o \
	module_nst_model.o \
	module_nst_parameters.o \
	module_nst_water_prop.o \
	moninp.o \
	moninp1.o \
	moninq.o \
	mstadb.o \
	mstadbtn.o \
	mstcnv.o \
	namelist_soilveg.o \
	ozphys.o \
	ozne_def.o \
	physcons.o \
	physparam.o \
	precpd.o \
	precpdp.o \
	progt2.o \
	progtm_module.o \
	radiation_aerosols.o \
	radiation_astronomy.o \
	radiation_clouds.o \
	radiation_gases.o \
	radiation_surface.o \
	radlw_datatb.o \
	radlw_main.o \
	radlw_param.o \
	radsw_datatb.o \
	radsw_main.o \
	radsw_param.o \
	rayleigh_damp.o \
	rayleigh_damp_mesopause.o \
	rascnvv2.o \
	sascnv.o \
	sascnvn.o \
	set_soilveg.o \
	sfc_diag.o \
	sfc_diff.o \
	sfc_drv.o \
	sfc_land.o \
	sfc_nst.o \
	sfc_ocean.o \
	sfc_sice.o \
	sflx.o \
	shalcnv.o \
	shalcv.o \
	shalcv_opr.o \
	tridi2t3.o

physics_nggps_src: $(OBJS)
	ar -ru ./../libphys.a $(OBJS)

# DEPENDENCIES:
calpreciptype.o		  : calpreciptype.f funcphys.o physcons.o
cnvc90.o		  : cnvc90.f machine.o
dcyc2.o			  : dcyc2.f machine.o physcons.o
dcyc2.pre.rad.o		  : dcyc2.pre.rad.f machine.o physcons.o
date_def.o		  : date_def.f machine.o
funcphys.o		  : funcphys.f machine.o physcons.o 
gbphys.o		  : gbphys.f machine.o physcons.o
get_prs.o		  : get_prs.f machine.o physcons.o tracer_const_h.o
gfs_phy_tracer_config.o   : gfs_phy_tracer_config.f machine.o
grrad.o			  : grrad.f machine.o funcphys.o physcons.o physparam.o radiation_aerosols.o radiation_astronomy.o \
				radiation_clouds.o radiation_gases.o radiation_surface.o radlw_main.o radlw_param.o \
				radsw_main.o radsw_param.o
gscond.o		  : gscond.f machine.o funcphys.o physcons.o
gscondp.o                 : gscondp.f machine.o funcphys.o physcons.o
gsmddrive.o		  : gsmddrive.f machine.o module_bfmicrophysics.o
gwdc.o			  : gwdc.f machine.o
gwdps.o			  : gwdps.f machine.o
lrgsclr.o		  : lrgsclr.f machine.o funcphys.o physcons.o
mersenne_twister.o        : mersenne_twister.f
module_bfmicrophysics.o	  : module_bfmicrophysics.f machine.o funcphys.o physcons.o
module_nst_model.o	  : module_nst_model.f module_nst_parameters.o module_nst_water_prop.o
module_nst_parameters.o	  : module_nst_parameters.f machine.o physcons.o
module_nst_water_prop.o	  : module_nst_water_prop.f machine.o module_nst_parameters.o
moninedmf.o               : moninedmf.f machine.o funcphys.o physcons.o
moninp.o		  : moninp.f machine.o physcons.o
moninp1.o		  : moninp1.f machine.o physcons.o
moninq.o		  : moninq.f machine.o funcphys.o physcons.o
mstadb.o		  : mstadb.f machine.o funcphys.o physcons.o
mstabdtn.o		  : mstabdtn.f machine.o funcphys.o physcons.o
mstcnv.o 		  : mstcnv.f machine.o funcphys.o physcons.o
namelist_soilveg.o	  : namelist_soilveg.f machine.o
ozne_def.o                : ozne_def.f machine.o
ozphys.o		  : ozphys.f machine.o physcons.o
physcons.o		  : physcons.f machine.o
physparam.o               : physparam.f machine.o
precpd.o		  : precpd.f machine.o funcphys.o physcons.o
precpdp.o		  : precpdp.f machine.o funcphys.o physcons.o
progt2.o		  : progt2.f machine.o progtm_module.o
progtm_module.o		  : progtm_module.f machine.o
radiation_aerosols.o	  : radiation_aerosols.f gfs_phy_tracer_config.o machine.o iounitdef.o physcons.o \
				physparam.o radlw_param.o radsw_param.o
radiation_astronomy.o	  : radiation_astronomy.f machine.o iounitdef.o physcons.o physparam.o
radiation_clouds.o	  : radiation_clouds.f machine.o iounitdef.o module_bfmicrophysics.o physcons.o physparam.o
radiation_gases.o	  : radiation_gases.f machine.o funcphys.o iounitdef.o ozne_def.o physcons.o physparam.o
radiation_surface.o	  : radiation_surface.f machine.o iounitdef.o physcons.o physparam.o
radlw_datatb.o		  : radlw_datatb.f machine.o physparam.o radlw_param.o
radlw_main.o		  : radlw_main.f machine.o mersenne_twister.o physparam.o radlw_datatb.o radlw_param.o
radlw_param.o		  : radlw_param.f machine.o physparam.o
radsw_datatb.o		  : radsw_datatb.f machine.o physparam.o radsw_param.o
radsw_main.o		  : radsw_main.f machine.o physcons.o physparam.o radsw_datatb.o radsw_param.o
radsw_param.o		  : radsw_param.f machine.o physparam.o
rascnvv2.o		  : rascnvv2.f machine.o funcphys.o
rayleigh_damp.o           : rayleigh_damp.f machine.o
rayleigh_damp_menopause.o : rayleigh_damp_menopause.f machine.o
sascnv.o		  : sascnv.f machine.o funcphys.o physcons.o
sascnvn.o		  : sascnvn.f machine.o funcphys.o physcons.o
set_soilveg.o		  : set_soilveg.f namelist_soilveg.o
sfc_diag.o		  : sfc_diag.f machine.o funcphys.o physcons.o
sfc_diff.o		  : sfc_diff.f machine.o funcphys.o physcons.o
sfc_drv.o		  : sfc_drv.f machine.o funcphys.o physcons.o
sfc_land.o		  : sfc_land.f machine.o funcphys.o physcons.o
sfc_nst.o		  : sfc_nst.f date_def.o machine.o funcphys.o module_nst_parameters.o module_nst_water_prop.o \
				physcons.o
sfc_ocean.o		  : sfc_ocean.f machine.o funcphys.o physcons.o
sfc_sice.o 		  : sfc_sice.f machine.o funcphys.o physcons.o
sflx.o			  : sflx.f machine.o physcons.o
shalcnv.o		  : shalcnv.f machine.o funcphys.o physcons.o
shalcv.o		  : shalcv.f machine.o physcons.o
shalcv_opr.o		  : shalcv_opr.f machine.o physcons.o
tracer_const_h.o	  : tracer_const_h.f machine.o physcons.o
tridi2t3.o		  : tridi2t3.f machine.o

cnvc90.o \
gscond.o \
gscondp.o \
get_prs.o \
gsmddrive.o \
gwdc.o \
gwdps.o \
lrgsclr.o \
mersenne_twister.o \
moninedmf.o \
moninp.o \
moninp1.o \
moninq.o \
mstadb.o \
mstadbtn.o \
mstcnv.o \
ozne_def.o \
ozphys.o \
precpd.o \
precpdp.o \
progt2.o \
progtm_module.o \
rayleigh_damp.o \
rayleigh_damp_mesopause.o \
sascnv.o \
sascnvn.o \
set_soilveg.o \
sfc_diag.o \
sfc_diff.o \
sfc_nst.o \
shalcnv.o \
shalcv.o \
shalcv_opr.o \
tracer_const_h.o \
tridi2t3.o:
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) $*.f > $*.f90
	$(FC) $(FFLAGS77_OPT) -c $*.f90 $(FCINCLUDES) -I../../../framework -I../../../operators -I..
else
	$(FC) $(CPPFLAGS) $(COREDEF) $(FFLAGS77_OPT) -c $*.f $(CPPINCLUDES) $(FCINCLUDES) -I../../../framework -I../../../operators -I..
endif

.f.o:
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../../../framework -I../../../operators -I..
else
	$(FC) $(CPPFLAGS) $(COREDEF) $(FFLAGS) -c $*.f $(CPPINCLUDES) $(FCINCLUDES) -I../../../framework -I../../../operators -I..
endif

.F.o:
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../../../framework -I../../../operators -I..
else
	$(FC) $(CPPFLAGS) $(COREDEF) $(FFLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES) -I../../../framework -I../../../operators -I..
endif

clean:
	$(RM) *.f90 *.o *.mod

