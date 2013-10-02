.SUFFIXES: .f .o

all: dummy physics_gfs_lib

dummy:
	echo "****** compile physics_gfs_lib ******"

OBJS = \
	iw3jdn.o \
	w3fs26.o \
	w3movdat.o \
	w3reddat.o

physics_gfs_lib: $(OBJS)
	ar -ru ./../libphys.a $(OBJS)

iw3jdn.o:
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) -DIWORDSIZE=4 -DRWORDSIZE=8 iw3jdn.f > iw3jdn.f90
	$(FC) $(FFLAGS77_OPT) -c iw3jdn.f90 $(FCINCLUDES) -I../../framework -I../../operators -I..

w3fs26.o:
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) -DIWORDSIZE=4 -DRWORDSIZE=8 w3fs26.f > w3fs26.f90
	$(FC) $(FFLAGS77_OPT) -c w3fs26.f90 $(FCINCLUDES) -I../../framework -I../../operators -I..

w3reddat.o:
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) -DIWORDSIZE=4 -DRWORDSIZE=8 w3reddat.f > w3reddat.f90
	$(FC) $(FFLAGS77_OPT) -c w3reddat.f90 $(FCINCLUDES) -I../../framework -I../../operators -I..

clean:
	$(RM) *.f90 *.o *.mod

.f.o:
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../../framework -I../../operators -I..

