.SUFFIXES: .f .o

all: dummy physics_gfs

dummy:
	echo "****** compile physics_gfs_src ******"

OBJS = \
	NameList_ESMFMod.o \
	Nst_Var_ESMFMod.o \
	Sfc_Var_ESMFMod.o \
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
	fix_fields.o \
	funcphys.o \
	gbphys.o \
	gcycle.o \
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
	nstio_module.o \
	ozinterp.o \
	ozphys.o \
	para_nstio_w.o \
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
	sfcio_module.o \
	sfcsub.o \
	sflx.o \
	shalcnv.o \
	shalcv.o \
	shalcv_opr.o \
	sig_io.o \
	tracer_const_h-new.o \
	tridi2t3.o \
	vert_def.o \
	wrtnst.o

physics_gfs: $(OBJS)
	ar -ru ./../libphys.a $(OBJS)

# DEPENDENCIES:
Nst_Var_ESMFMod.o	: machine.o
Sfc_Var_EMSFMod.o	: machine.o
atm.comm.o		: machine.o layout1.o physcons.o mpi_def.o namelist_def.o
calpreciptype.o		: funcphys.o physcons.o
cmp.comm.o		: machine.o
cnvc90.o		: machine.o
compns.o		: layout_grid_tracers.o mpi_def.o namelist_def.o
coordinate_def.o	: machine.o
d3d_def.o		: machine.o
dcyc2.o			: machine.o physcons.o
dcyc2.pre.rad.o		: machine.o physcons.o
date_def.o		: machine.o
fix_fields.o		: machine.o funcphys.o layout1.o progtm_module.o namelist_def.o resol_def.o \
				Nst_Var_ESMFMod.o Sfc_Var_ESMFMod.o 
funcphys.o		: machine.o physcons.o 
gbphys.o		: machine.o physcons.o
gcycle.o		: machine.o physcons.o resol_def.o Sfc_Var_ESMFMod.o
get_prs.o		: machine.o physcons.o tracer_const_h-new.o
gfsmisc_def.o           : machine.o
gg_def.o		: machine.o
grrad.o			: machine.o funcphys.o physcons.o radiation_aerosols.o radiation_astronomy.o radiation_clouds.o \
				radiation_gases.o radiation_surface.o radlw_main.o radlw_param.o radsw_main.o radsw_param.o
gscond.o		: machine.o funcphys.o physcons.o
gsmddrive.o		: machine.o module_bfmicrophysics.o
gwdc.o			: machine.o
gwdps.o			: machine.o
lrgsclr.o		: machine.o funcphys.o physcons.o
layout_grid_tracers.o	: machine.o
matrix_sig_def.o	: machine.o
module_bfmicrophysics.o	: machine.o funcphys.o physcons.o
module_nst_model.o	: module_nst_parameters.o module_nst_water_prop.o
module_nst_parameters.o	: machine.o physcons.o
module_nst_water_prop.o	: machine.o module_nst_parameters.o
moninp.o		: machine.o physcons.o
moninp1.o		: machine.o physcons.o
moninq.o		: machine.o funcphys.o physcons.o
moninq1.o		: machine.o funcphys.o physcons.o
mpi_def.o		: machine.o
mpi_quit.o		: mpi_def.o
mstadb.o		: machine.o funcphys.o physcons.o
mstabdtn.o		: machine.o funcphys.o physcons.o
mstcnv.o 		: machine.o funcphys.o physcons.o
namelist_def.o		: machine.o
namelist_soilveg.o	: machine.o
ozinterp.o		: machine.o gg_def.o layout1.o resol_def.o
ozphys.o		: machine.o physcons.o
para_nstio_w.o		: resol_def.o layout1.o namelist_def.o nstio_module.o Nst_Var_ESMFMod.o
physcons.o		: machine.o
precpd.o		: machine.o funcphys.o physcons.o
progt2.o		: machine.o progtm_module.o
progtm_module.o		: machine.o
radiation_aerosols.o	: machine.o iounitdef.o physcons.o radlw_param.o radsw_param.o
radiation_astronomy.o	: machine.o iounitdef.o physcons.o
radiation_clouds.o	: machine.o iounitdef.o module_bfmicrophysics.o physcons.o
radiation_gases.o	: machine.o funcphys.o iounitdef.o physcons.o resol_def.o
radiation_surface.o	: machine.o iounitdef.o physcons.o
radlw_datatb.o		: machine.o radlw_param.o
radlw_main.o		: machine.o radlw_datatb.o radlw_param.o
radlw_param.o		: machine.o
radsw_datatb.o		: machine.o radsw_param.o
radsw_main.o		: machine.o physcons.o radsw_datatb.o radsw_param.o
radsw_param.o		: machine.o
rascnvv2.o		: machine.o funcphys.o
resol_def.o		: machine.o
sascnv.o		: machine.o funcphys.o physcons.o
sascnvn.o		: machine.o funcphys.o physcons.o
set_soilveg.o		: namelist_soilveg.o
sfc_diag.o		: machine.o funcphys.o physcons.o
sfc_diff.o		: machine.o funcphys.o physcons.o
sfc_drv.o		: machine.o funcphys.o physcons.o
sfc_land.o		: machine.o funcphys.o physcons.o
sfc_nst.o		: date_def.o machine.o funcphys.o module_nst_parameters.o module_nst_water_prop.o \
				physcons.o
sfc_ocean.o		: machine.o funcphys.o physcons.o
sfc_sice.o 		: machine.o funcphys.o physcons.o
sfcsub.o		: machine.o
sflx.o			: machine.o physcons.o
shalcnv.o		: machine.o funcphys.o physcons.o
shalcv.o		: machine.o physcons.o
shalcv_opr.o		: machine.o physcons.o
sig_io.o		: resol_def.o
tracer_const_h-new.o	: machine.o physcons.o
tridi2t3.o		: machine.o
vert_def.o		: machine.o
wrtnst.o		: machine.o layout1.o namelist_def.o resol_def.o sig_io.o

clean:
	$(RM) *.f90 *.o *.mod

atm.comm.o:
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) cmp.comm.f > atm.comm.f90
	$(FC) $(FFLAGS77_OPT) -c atm.comm.f90 $(FCINCLUDES) -I../../../framework -I../../../operators -I..

cnvc90.o:
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) cnvc90.f > cnvc90.f90
	$(FC) $(FFLAGS77_OPT) -c cnvc90.f90 $(FCINCLUDES) -I../../../framework -I../../../operators -I..

cmp.comm.o:
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) -DIWORDSIZE=4 -DRWORDSIZE=8 cmp.comm.f > cmp.comm.f90
	$(FC) $(FFLAGS77_OPT) -c cmp.comm.f90 $(FCINCLUDES) -I../../../framework -I../../../operators -I..

compns.o:
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) -DIWORDSIZE=4 -DRWORDSIZE=8 compns.f > compns.f90
	$(FC) $(FFLAGS77_OPT) -c compns.f90 $(FCINCLUDES) -I../../../framework -I../../../operators -I..

coordinate_def.o:
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) -DIWORDSIZE=4 -DRWORDSIZE=8 coordinate_def.f > coordinate_def.f90
	$(FC) $(FFLAGS77_OPT) -c coordinate_def.f90 $(FCINCLUDES) -I../../../framework -I../../../operators -I..

d3d_def.o:
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) -DIWORDSIZE=4 -DRWORDSIZE=8 d3d_def.f > d3d_def.f90
	$(FC) $(FFLAGS77_OPT) -c d3d_def.f90 $(FCINCLUDES) -I../../../framework -I../../../operators -I..

fix_fields.o:
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) -DIWORDSIZE=4 -DRWORDSIZE=8 fix_fields.f > fix_fields.f90
	$(FC) $(FFLAGS77_OPT) -c fix_fields.f90 $(FCINCLUDES) -I../../../framework -I../../../operators -I..

gcycle.o:
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) -DIWORDSIZE=4 -DRWORDSIZE=8 gcycle.f > gcycle.f90
	$(FC) $(FFLAGS77_OPT) -c gcycle.f90 $(FCINCLUDES) -I../../../framework -I../../../operators -I..

gscond.o:
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) -DIWORDSIZE=4 -DRWORDSIZE=8 gscond.f > gscond.f90
	$(FC) $(FFLAGS77_OPT) -c gscond.f90 $(FCINCLUDES) -I../../../framework -I../../../operators -I..

get_prs.o:
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) -DIWORDSIZE=4 -DRWORDSIZE=8 get_prs.f > get_prs.f90
	$(FC) $(FFLAGS77_OPT) -c get_prs.f90 $(FCINCLUDES) -I../../../framework -I../../../operators -I..

gfsmisc_def.o:
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) -DIWORDSIZE=4 -DRWORDSIZE=8 gfsmisc_def.f > gfsmisc_def.f90
	$(FC) $(FFLAGS77_OPT) -c gfsmisc_def.f90 $(FCINCLUDES) -I../../../framework -I../../../operators -I..

gg_def.o:
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) -DIWORDSIZE=4 -DRWORDSIZE=8 gg_def.f > gg_def.f90
	$(FC) $(FFLAGS77_OPT) -c gg_def.f90 $(FCINCLUDES) -I../../../framework -I../../../operators -I..

gwdc.o:
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) -DIWORDSIZE=4 -DRWORDSIZE=8 gwdc.f > gwdc.f90
	$(FC) $(FFLAGS77_OPT) -c gwdc.f90 $(FCINCLUDES) -I../../../framework -I../../../operators -I..

gwdps.o:
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) -DIWORDSIZE=4 -DRWORDSIZE=8 gwdps.f > gwdps.f90
	$(FC) $(FFLAGS77_OPT) -c gwdps.f90 $(FCINCLUDES) -I../../../framework -I../../../operators -I..

layout1.o:
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) -DIWORDSIZE=4 -DRWORDSIZE=8 layout1.f > layout1.f90
	$(FC) $(FFLAGS77_OPT) -c layout1.f90 $(FCINCLUDES) -I../../../framework -I../../../operators -I..

layout_grid_tracers.o:
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) -DIWORDSIZE=4 -DRWORDSIZE=8 layout_grid_tracers.f > layout_grid_tracers.f90
	$(FC) $(FFLAGS77_OPT) -c layout_grid_tracers.f90 $(FCINCLUDES) -I../../../framework -I../../../operators -I..

lrgsclr.o:
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) -DIWORDSIZE=4 -DRWORDSIZE=8 lrgsclr.f > lrgsclr.f90
	$(FC) $(FFLAGS77_OPT) -c lrgsclr.f90 $(FCINCLUDES) -I../../../framework -I../../../operators -I..

machine.o:
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) -DIWORDSIZE=4 -DRWORDSIZE=8 machine.f > machine.f90
	$(FC) $(FFLAGS77_OPT) -c machine.f90 $(FCINCLUDES) -I../../../framework -I../../../operators -I..

matrix_sig_def.o:
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) -DIWORDSIZE=4 -DRWORDSIZE=8 matrix_sig_def.f > matrix_sig_def.f90
	$(FC) $(FFLAGS77_OPT) -c matrix_sig_def.f90 $(FCINCLUDES) -I../../../framework -I../../../operators -I..

mersenne_twister.o:
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) -DIWORDSIZE=4 -DRWORDSIZE=8 mersenne_twister.f > mersenne_twister.f90
	$(FC) $(FFLAGS77_OPT) -c mersenne_twister.f90 $(FCINCLUDES) -I../../../framework -I../../../operators -I..

moninp.o:
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) -DIWORDSIZE=4 -DRWORDSIZE=8 moninp.f > moninp.f90
	$(FC) $(FFLAGS77_OPT) -c moninp.f90 $(FCINCLUDES) -I../../../framework -I../../../operators -I..

moninp1.o:
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) -DIWORDSIZE=4 -DRWORDSIZE=8 moninp1.f > moninp1.f90
	$(FC) $(FFLAGS77_OPT) -c moninp1.f90 $(FCINCLUDES) -I../../../framework -I../../../operators -I..

moninq.o:
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) -DIWORDSIZE=4 -DRWORDSIZE=8 moninq.f > moninq.f90
	$(FC) $(FFLAGS77_OPT) -c moninq.f90 $(FCINCLUDES) -I../../../framework -I../../../operators -I..

moninq1.o:
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) -DIWORDSIZE=4 -DRWORDSIZE=8 moninq1.f > moninq1.f90
	$(FC) $(FFLAGS77_OPT) -c moninq1.f90 $(FCINCLUDES) -I../../../framework -I../../../operators -I..

mpi_def.o:
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) -DIWORDSIZE=4 -DRWORDSIZE=8 mpi_def.f > mpi_def.f90
	$(FC) $(FFLAGS77_OPT) -c mpi_def.f90 $(FCINCLUDES) -I../../../framework -I../../../operators -I..

mpi_quit.o:
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) -DIWORDSIZE=4 -DRWORDSIZE=8 mpi_quit.f > mpi_quit.f90
	$(FC) $(FFLAGS77_OPT) -c mpi_quit.f90 $(FCINCLUDES) -I../../../framework -I../../../operators -I..

mstadb.o:
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) -DIWORDSIZE=4 -DRWORDSIZE=8 mstadb.f > mstadb.f90
	$(FC) $(FFLAGS77_OPT) -c mstadb.f90 $(FCINCLUDES) -I../../../framework -I../../../operators -I..

mstadbtn.o:
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) -DIWORDSIZE=4 -DRWORDSIZE=8 mstadbtn.f > mstadbtn.f90
	$(FC) $(FFLAGS77_OPT) -c mstadbtn.f90 $(FCINCLUDES) -I../../../framework -I../../../operators -I..

mstcnv.o:
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) -DIWORDSIZE=4 -DRWORDSIZE=8 mstcnv.f > mstcnv.f90
	$(FC) $(FFLAGS77_OPT) -c mstcnv.f90 $(FCINCLUDES) -I../../../framework -I../../../operators -I..

namelist_def.o:
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) -DIWORDSIZE=4 -DRWORDSIZE=8 namelist_def.f > namelist_def.f90
	$(FC) $(FFLAGS77_OPT) -c namelist_def.f90 $(FCINCLUDES) -I../../../framework -I../../../operators -I..

ozinterp.o:
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) -DIWORDSIZE=4 -DRWORDSIZE=8 ozinterp.f > ozinterp.f90
	$(FC) $(FFLAGS77_OPT) -c ozinterp.f90 $(FCINCLUDES) -I../../../framework -I../../../operators -I..

ozphys.o:
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) -DIWORDSIZE=4 -DRWORDSIZE=8 ozphys.f > ozphys.f90
	$(FC) $(FFLAGS77_OPT) -c ozphys.f90 $(FCINCLUDES) -I../../../framework -I../../../operators -I..

para_nstio_w.o:
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) -DIWORDSIZE=4 -DRWORDSIZE=8 para_nstio_w.f > para_nstio_w.f90
	$(FC) $(FFLAGS77_OPT) -c para_nstio_w.f90 $(FCINCLUDES) -I../../../framework -I../../../operators -I..

precpd.o:
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) -DIWORDSIZE=4 -DRWORDSIZE=8 precpd.f > precpd.f90
	$(FC) $(FFLAGS77_OPT) -c precpd.f90 $(FCINCLUDES) -I../../../framework -I../../../operators -I..

progt2.o:
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) -DIWORDSIZE=4 -DRWORDSIZE=8 progt2.f > progt2.f90
	$(FC) $(FFLAGS77_OPT) -c progt2.f90 $(FCINCLUDES) -I../../../framework -I../../../operators -I..

progtm_module.o:
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) -DIWORDSIZE=4 -DRWORDSIZE=8 progtm_module.f > progtm_module.f90
	$(FC) $(FFLAGS77_OPT) -c progtm_module.f90 $(FCINCLUDES) -I../../../framework -I../../../operators -I..

resol_def.o:
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) -DIWORDSIZE=4 -DRWORDSIZE=8 resol_def.f > resol_def.f90
	$(FC) $(FFLAGS77_OPT) -c resol_def.f90 $(FCINCLUDES) -I../../../framework -I../../../operators -I..

sascnv.o:
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) -DIWORDSIZE=4 -DRWORDSIZE=8 sascnv.f > sascnv.f90
	$(FC) $(FFLAGS77_OPT) -c sascnv.f90 $(FCINCLUDES) -I../../../framework -I../../../operators -I..

sascnvn.o:
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) -DIWORDSIZE=4 -DRWORDSIZE=8 sascnvn.f > sascnvn.f90
	$(FC) $(FFLAGS77_OPT) -c sascnvn.f90 $(FCINCLUDES) -I../../../framework -I../../../operators -I..

set_soilveg.o:
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) -DIWORDSIZE=4 -DRWORDSIZE=8 set_soilveg.f > set_soilveg.f90
	$(FC) $(FFLAGS77_OPT) -c set_soilveg.f90 $(FCINCLUDES) -I../../../framework -I../../../operators -I..

sfc_diag.o:
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) -DIWORDSIZE=4 -DRWORDSIZE=8 sfc_diag.f > sfc_diag.f90
	$(FC) $(FFLAGS77_OPT) -c sfc_diag.f90 $(FCINCLUDES) -I../../../framework -I../../../operators -I..

sfc_diff.o:
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) -DIWORDSIZE=4 -DRWORDSIZE=8 sfc_diff.f > sfc_diff.f90
	$(FC) $(FFLAGS77_OPT) -c sfc_diff.f90 $(FCINCLUDES) -I../../../framework -I../../../operators -I..

sfc_nst.o:
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) -DIWORDSIZE=4 -DRWORDSIZE=8 sfc_nst.f > sfc_nst.f90
	$(FC) $(FFLAGS77_OPT) -c sfc_nst.f90 $(FCINCLUDES) -I../../../framework -I../../../operators -I..

sfcsub.o:
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) -DIWORDSIZE=4 -DRWORDSIZE=8 sfcsub.f > sfcsub.f90
	$(FC) $(FFLAGS77_OPT) -c sfcsub.f90 $(FCINCLUDES) -I../../../framework -I../../../operators -I..

shalcnv.o:
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) -DIWORDSIZE=4 -DRWORDSIZE=8 shalcnv.f > shalcnv.f90
	$(FC) $(FFLAGS77_OPT) -c shalcnv.f90 $(FCINCLUDES) -I../../../framework -I../../../operators -I..

shalcv.o:
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) -DIWORDSIZE=4 -DRWORDSIZE=8 shalcv.f > shalcv.f90
	$(FC) $(FFLAGS77_OPT) -c shalcv.f90 $(FCINCLUDES) -I../../../framework -I../../../operators -I..

shalcv_opr.o:
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) -DIWORDSIZE=4 -DRWORDSIZE=8 shalcv_opr.f > shalcv_opr.f90
	$(FC) $(FFLAGS77_OPT) -c shalcv_opr.f90 $(FCINCLUDES) -I../../../framework -I../../../operators -I..

sig_io.o:
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) -DIWORDSIZE=4 -DRWORDSIZE=8 sig_io.f > sig_io.f90
	$(FC) $(FFLAGS77_OPT) -c sig_io.f90 $(FCINCLUDES) -I../../../framework -I../../../operators -I..

tridi2t3.o:
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) -DIWORDSIZE=4 -DRWORDSIZE=8 tridi2t3.f > tridi2t3.f90
	$(FC) $(FFLAGS77_OPT) -c tridi2t3.f90 $(FCINCLUDES) -I../../../framework -I../../../operators -I..

vert_def.o:
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) -DIWORDSIZE=4 -DRWORDSIZE=8 vert_def.f > vert_def.f90
	$(FC) $(FFLAGS77_OPT) -c vert_def.f90 $(FCINCLUDES) -I../../../framework -I../../../operators -I..

wrtnst.o:
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) -DIWORDSIZE=4 -DRWORDSIZE=8 wrtnst.f > wrtnst.f90
	$(FC) $(FFLAGS77_OPT) -c wrtnst.f90 $(FCINCLUDES) -I../../../framework -I../../../operators -I..

.f.o:
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) -DIWORDSIZE=4 -DRWORDSIZE=8 $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I./../physics_gfs_lib -I../../../framework -I../../../operators -I..

