
FC = ifort
FFLAGS = -O2 -132 -traceback
LDFLAGS = -mkl

src = thermostat.o const.o ctrl.o main.o initialize.o simulation.o NHC_VV.o MD_VV.o
pes = pes_pbe_Ru.o pes_pbe_int.o pes_tbody_int.o pes_eann.o \
      pes_interface3.o \
      energy1_int.o pes_ococ_v4.2.o \
      energy1_co.o

vmd.x:	$(pes) $(src)
	$(FC) $(FFLAGS) $+ -o $@ $(LDFLAGS)

%.o:	%.f90
	$(FC) $(FFLAGS) -c $< -o $@
%.o:	%.F90
	$(FC) $(FFLAGS) -c $< -o $@
%.o:	%.f
	$(FC) $(FFLAGS) -c $< -o $@
%.o:	%.F
	$(FC) $(FFLAGS) -c $< -o $@

