
FC = ifort
FFLAGS = -O2 -132
LDFLAGS = -mkl

src = thermostat.o const.o ctrl.o main.o initialize.o simulation.o NHC_VV.o MD_VV.o
#pes = pes_interface.f90 pes_interface2.f90 pes_ococ_v4.2.F90 pes_copt6d_vdw.f90 energy1_int.F 
pes = pes_hd.f90 pes_ruFixco.f90 pes_interfaceru.f90 

vmd.x:	$(src) $(pes)
	$(FC) $(FFLAGS) $+ -o $@ $(LDFLAGS)

%.o:	%.f90
	$(FC) $(FFLAGS) -c $< -o $@
%.o:	%.F90
	$(FC) $(FFLAGS) -c $< -o $@
%.o:	%.f
	$(FC) $(FFLAGS) -c $< -o $@
%.o:	%.F
	$(FC) $(FFLAGS) -c $< -o $@
clean:
	rm *.mod *.o vmd.x

