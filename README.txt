

The MD Code for CO/Ru(0001) System

MD_Code: 
Umbrellasampling:
the code of umbrella sampling MD simulation for CO desorption using EANN potential or DPA2 fine-tuned potential, MD simulation for EANN (2x2x5 and 4x4x5) potential using Fortran and DPA2 fine-tuned potential (20x20x5) using python.
UnbiasedMD:
the code of  MD simulation for CO diffusion on Ru(0001) surface using EANN potential.

ML_Models: the PES interface and PES file for each CO/Ru(0001)  System.
2x2x5_EANN contains 6D PES and 42D pes
the 6D pes contains 2 files, "nnfit_3.14.txt" and "pes_ruFIXco.f90"
the HD pes contains 2 files and 1 folder, "pes_hd.f90", "pes_hd.inc1.f90", and "pes8.dat"

4x4x5_EANN is a high-dimensional PES for CO on 4x4x5 Ru(0001) surface constructed by many body expansion strategy, which contains 8 files and 1 foler, pes_eann.f90, pes_interface3.f90, pes_ococ_v4.2.F90, pes_ococ_w87.txt, pes_pbe_int.f90, pes_pbe_Ru.f90, pes_tbody_int.f90, and 4x4x5_pes_coru

20x20x5_DPA2 is a DPA2 fine-tuned PES for CO/Ru(0001) for performing large-cell MD simulation, which contains 1 file, dpa2_m_OC20M247f.pth

Example:
input-445(folder) : inputs for umbrella sampling MD simulations for 0.500 ML on 4x4x5 Ru(0001) surface
input-20205(folder) : inputs for umbrella sampling MD simulations for 0.500 ML on 20x20x5 Ru(0001) surface

Umbrella Sampling:
The total code of each size umbrella sampling for CO/Ru(0001) system, complie using "make" command when using  EANN PES (2x2x5 and 4x4x5) and using python dirofmain/main.py for DPA2 fine-tuned PES.
US-225:  umbrella sampling MD simulation code for CO on 2x2x5 Ru(0001) surface, INPUT is the sample of 6DPES and INPUT.2layer for 42D PES.
US-445: umbrella sampling MD simulation code for CO on 4x4x5 Ru(0001) surface, INPUT is the sample of two COs on two layer relaxed surface. The "Example" folder provided more explicit input files for each window.
US-DPA2: umbrella sampling MD simulation code for CO on 20x20x5 Ru(0001) surface, INPUT is the sample for a window of 0.500 ML coverage on 20x20 Ru(0001) surface. The version of deepmd-kit is 3.1.1 with CUDA 12.9.

UnbiasedMD:
The total code of each size umbrella sampling for CO/Ru(0001) system, using "make" command to complie
2x2x5_EANN:the code of  MD simulation for CO diffusion on 2x2x5 Ru(0001) surface using EANN potential. INPUT is the sample of one CO on two layer relaxed surface.
4x4x5_EANN:the code of  MD simulation for CO diffusion on 4x4x5 Ru(0001) surface using EANN potential. INPUT is the sample of three COs on two layer relaxed surface.

Dataset:
The training data set for 4x4x5 EANN potential and dpa2 fine-tuned  PES
4x4x5_surf(folder) : The training data set for 4x4x5 Ru(0001) surface PES
4x4x5_rucoint(folder) : The training data set for CO-Ru(0001) adsorbate-surface interaction PES
4x4x5_surf(folder) : The training data set for CO-CO-Ru(0001) threeboy interaction PES
dpa2-finetune(folder) : The training data set fordpa2 fine-tuned  PES
