.SUFFIXES: .f .o

all: dummy physics_gfs_src

dummy:
	echo "****** compile physics_gfs_src ******"

OBJS = \
	atm.comm.o \
	calpreciptype.o \
	cmp.comm.o \
	compns.o \
	coordinate_def.o \
	cnvc90.o \
	d3d_def.o \
	date_def.o \
	dcyc2.o \
	dcyc2.pre.rad.o \
	funcphys.o \
	gbphys.o \
	get_prs.o \
	gfsmisc_def.o \
	gg_def.o \
	grrad.o \
	gscond.o \
	gwdc.o \
	gwdps.o \
	gsmddrive.o \
        iounitdef.o \
	layout1.o \
	lrgsclr.o \
	layout_grid_tracers.o \
	machine.o  \
	matrix_sig_def.o \
	mersenne_twister.o \
	module_bfmicrophysics.o \
	module_nst_model.o \
	module_nst_parameters.o \
	module_nst_water_prop.o \
	moninp.o \
	moninp1.o \
	moninq.o \
	moninq1.o \
	mpi_def.o \
	mpi_quit.o \
	mstadb.o \
	mstadbtn.o \
	mstcnv.o \
	namelist_def.o \
	namelist_soilveg.o \
	ozphys.o \
	precpd.o \
	progt2.o \
	progtm_module.o \
	physcons.o \
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
	rascnvv2.o \
	resol_def.o \
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
	tracer_const_h-new.o \
	tridi2t3.o \
	vert_def.o

physics_gfs_src: $(OBJS)
	ar -ru ./../libphys.a $(OBJS)

# DEPENDENCIES:
atm.comm.o		: atm.comm.f machine.o cmp.comm.o layout1.o physcons.o mpi_def.o namelist_def.o
calpreciptype.o		: calpreciptype.f funcphys.o physcons.o
cmp.comm.o		: cmp.comm.f machine.o
cnvc90.o		: cnvc90.f machine.o
compns.o		: compns.f layout_grid_tracers.o mpi_def.o namelist_def.o
coordinate_def.o	: coordinate_def.f machine.o
d3d_def.o		: d3d_def.f machine.o
dcyc2.o			: dcyc2.f machine.o physcons.o
dcyc2.pre.rad.o		: dcyc2.pre.rad.f machine.o physcons.o
date_def.o		: date_def.f machine.o
funcphys.o		: funcphys.f machine.o physcons.o 
gbphys.o		: gbphys.f machine.o physcons.o
get_prs.o		: get_prs.f machine.o physcons.o tracer_const_h-new.o
gfsmisc_def.o           : gfsmisc_def.f machine.o
gg_def.o		: gg_def.f machine.o
grrad.o			: grrad.f machine.o funcphys.o physcons.o radiation_aerosols.o radiation_astronomy.o \
				radiation_clouds.o radiation_gases.o radiation_surface.o radlw_main.o radlw_param.o \
				radsw_main.o radsw_param.o
gscond.o		: gscond.f machine.o funcphys.o physcons.o
gsmddrive.o		: gsmddrive.f machine.o module_bfmicrophysics.o
gwdc.o			: gwdc.f machine.o
gwdps.o			: gwdps.f machine.o
lrgsclr.o		: lrgsclr.f machine.o funcphys.o physcons.o
layout_grid_tracers.o	: layout_grid_tracers.f machine.o
matrix_sig_def.o	: matrix_sig_def.f machine.o
mersenne_twister.o      : mersenne_twister.f
module_bfmicrophysics.o	: module_bfmicrophysics.f machine.o funcphys.o physcons.o
module_nst_model.o	: module_nst_model.f module_nst_parameters.o module_nst_water_prop.o
module_nst_parameters.o	: module_nst_parameters.f machine.o physcons.o
module_nst_water_prop.o	: module_nst_water_prop.f machine.o module_nst_parameters.o
moninp.o		: moninp.f machine.o physcons.o
moninp1.o		: moninp1.f machine.o physcons.o
moninq.o		: moninq.f machine.o funcphys.o physcons.o
moninq1.o		: moninq1.f machine.o funcphys.o physcons.o
mpi_def.o		: mpi_def.f machine.o
mpi_quit.o		: mpi_quit.f mpi_def.o
mstadb.o		: mstadb.f machine.o funcphys.o physcons.o
mstabdtn.o		: mstabdtn.f machine.o funcphys.o physcons.o
mstcnv.o 		: mstcnv.f machine.o funcphys.o physcons.o
namelist_def.o		: namelist_def.f machine.o
namelist_soilveg.o	: namelist_soilveg.f machine.o
ozphys.o		: ozphys.f machine.o physcons.o
physcons.o		: physcons.f machine.o
precpd.o		: precpd.f machine.o funcphys.o physcons.o
progt2.o		: progt2.f machine.o progtm_module.o
progtm_module.o		: progtm_module.f machine.o
radiation_aerosols.o	: radiation_aerosols.f machine.o iounitdef.o physcons.o radlw_param.o radsw_param.o
radiation_astronomy.o	: radiation_astronomy.f machine.o iounitdef.o physcons.o
radiation_clouds.o	: radiation_clouds.f machine.o iounitdef.o module_bfmicrophysics.o physcons.o
radiation_gases.o	: radiation_gases.f machine.o funcphys.o iounitdef.o physcons.o resol_def.o
radiation_surface.o	: radiation_surface.f machine.o iounitdef.o physcons.o
radlw_datatb.o		: radlw_datatb.f machine.o radlw_param.o
radlw_main.o		: radlw_main.f machine.o radlw_datatb.o radlw_param.o
radlw_param.o		: radlw_param.f machine.o
radsw_datatb.o		: radsw_datatb.f machine.o radsw_param.o
radsw_main.o		: radsw_main.f machine.o physcons.o radsw_datatb.o radsw_param.o
radsw_param.o		: radsw_param.f machine.o
rascnvv2.o		: rascnvv2.f machine.o funcphys.o
resol_def.o		: resol_def.f machine.o
sascnv.o		: sascnv.f machine.o funcphys.o physcons.o
sascnvn.o		: sascnvn.f machine.o funcphys.o physcons.o
set_soilveg.o		: set_soilveg.f namelist_soilveg.o
sfc_diag.o		: sfc_diag.f machine.o funcphys.o physcons.o
sfc_diff.o		: sfc_diff.f machine.o funcphys.o physcons.o
sfc_drv.o		: sfc_drv.f machine.o funcphys.o physcons.o
sfc_land.o		: sfc_land.f machine.o funcphys.o physcons.o
sfc_nst.o		: sfc_nst.f date_def.o machine.o funcphys.o module_nst_parameters.o module_nst_water_prop.o \
				physcons.o
sfc_ocean.o		: sfc_ocean.f machine.o funcphys.o physcons.o
sfc_sice.o 		: sfc_sice.f machine.o funcphys.o physcons.o
sflx.o			: sflx.f machine.o physcons.o
shalcnv.o		: shalcnv.f machine.o funcphys.o physcons.o
shalcv.o		: shalcv.f machine.o physcons.o
shalcv_opr.o		: shalcv_opr.f machine.o physcons.o
tracer_const_h-new.o	: tracer_const_h-new.f machine.o physcons.o
tridi2t3.o		: tridi2t3.f machine.o
vert_def.o		: vert_def.f machine.o

atm.comm.o \
cnvc90.o \
cmp.comm.o \
compns.o \
coordinate_def.o \
d3d_def.o \
gscond.o \
get_prs.o \
gfsmisc_def.o \
gg_def.o \
gwdc.o \
gwdps.o \
layout1.o \
layout_grid_tracers.o \
lrgsclr.o \
matrix_sig_def.o \
mersenne_twister.o \
moninp.o \
moninp1.o \
moninq.o \
moninq1.o \
mpi_def.o \
mpi_quit.o \
mstadb.o \
mstadbtn.o \
mstcnv.o \
namelist_def.o \
ozphys.o \
precpd.o \
progt2.o \
progtm_module.o \
resol_def.o \
sascnv.o \
sascnvn.o \
set_soilveg.o \
sfc_diag.o \
sfc_diff.o \
sfc_nst.o \
shalcnv.o \
shalcv.o \
shalcv_opr.o \
tracer_const_h-new.o \
tridi2t3.o \
vert_def.o:
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) $*.f > $*.f90
	$(FC) $(FFLAGS77_OPT) -c $*.f90 $(FCINCLUDES) -I./../physics_gfs_lib -I../../../framework -I../../../operators -I..
else
	$(FC) $(CPPFLAGS) $(COREDEF) $(FFLAGS77_OPT) -c $*.f $(CPPINCLUDES) $(FCINCLUDES) -I./../physics_gfs_lib -I../../../framework -I../../../operators -I..
endif

.f.o:
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I./../physics_gfs_lib -I../../../framework -I../../../operators -I..
else
	$(FC) $(CPPFLAGS) $(COREDEF) $(FFLAGS) -c $*.f $(CPPINCLUDES) $(FCINCLUDES) -I./../physics_gfs_lib -I../../../framework -I../../../operators -I..
endif

clean:
	$(RM) *.f90 *.o *.mod

