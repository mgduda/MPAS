.SUFFIXES: .f .o

all: dummy physics_nggps_lib

dummy:
	echo "****** compile physics_nggps_lib ******"

OBJS = \
	iw3jdn.o \
	w3fs26.o \
	w3kind.o \
	w3movdat.o \
	w3reddat.o

physics_nggps_lib: $(OBJS)
	ar -ru ./../libphys.a $(OBJS)

ifeq "$(GEN_F90)" "true"

iw3jdn.o:
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) iw3jdn.f > iw3jdn.f90
	$(FC) $(FFLAGS77) -c iw3jdn.f90 $(FCINCLUDES) -I../../framework -I../../operators -I..
w3fs26.o:
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) w3fs26.f > w3fs26.f90
	$(FC) $(FFLAGS77) -c w3fs26.f90 $(FCINCLUDES) -I../../framework -I../../operators -I..
w3reddat.o:
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) w3reddat.f > w3reddat.f90
	$(FC) $(FFLAGS77) -c w3reddat.f90 $(FCINCLUDES) -I../../framework -I../../operators -I..
.f.o:
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../../framework -I../../operators -I..

else

iw3jdn.o:
	$(FC) $(CPPFLAGS) $(COREDEF) $(FFLAGS77) -c iw3jdn.f $(CPPINCLUDES) $(FCINCLUDES) -I../../framework -I../../operators -I..
w3fs26.o:
	$(FC) $(CPPFLAGS) $(COREDEF) $(FFLAGS77) -c w3fs26.f $(CPPINCLUDES) $(FCINCLUDES) -I../../framework -I../../operators -I..
w3reddat.o:
	$(FC) $(CPPFLAGS) $(COREDEF) $(FFLAGS77) -c w3reddat.f $(CPPINCLUDES) $(FCINCLUDES) -I../../framework -I../../operators -I..
.f.o:
	$(FC) $(CPPFLAGS) $(COREDEF) $(FFLAGS) -c $*.f $(CPPINCLUDES) $(FCINCLUDES) -I../../../framework -I../../../operators -I..

endif

clean:
	$(RM) *.f90 *.o *.mod

