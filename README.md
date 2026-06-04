# Source Code and Datasets for "Coverage-Dependent Free-Energy Bottlenecks in CO Adsorption and Desorption on Ru(0001)"

This repository contains the simulation source codes, machine-learning potential energy surfaces (ML-PES), example input profiles, data post-processing toolkits, and raw training datasets associated with our study on the coverage-dependent adsorption and desorption dynamics of CO on the Ru(0001) surface.

---

## Authors & Affiliations

* **Jin-En He** *State Key Laboratory of Structural Chemistry, Fujian Institute of Research on the Structure of Matter, Chinese Academy of Sciences, Fuzhou, 350002, China*
* **Mingjun Yang** *Shenzhen Jingtai Technology Co., Ltd. (XtalPi), Fubao Community, Shenzhen 518045, China*
* **Zhe-Ning Chen** *(Corresponding Author: znchen@fjirsm.ac.cn)* *State Key Laboratory of Structural Chemistry, Fujian Institute of Research on the Structure of Matter, Chinese Academy of Sciences, Fuzhou, 350002, China* *Fujian Key Laboratory of Theoretical and Computational Chemistry, Xiamen 361005, China*
* **Tong Zhu** *(Corresponding Author: tzhu@lps.ecnu.edu.cn)* *Shanghai Engineering Research Center of Molecular Therapeutics and New Drug Development, School of Chemistry and Molecular Engineering, East China Normal University, Shanghai 200062, China* *Shanghai Innovation Institute, Shanghai 200003, China*
* **Jun Chen** *(Corresponding Author: chenjun@fjirsm.ac.cn)* *State Key Laboratory of Structural Chemistry, Fujian Institute of Research on the Structure of Matter, Chinese Academy of Sciences, Fuzhou, 350002, China* *Fujian Key Laboratory of Theoretical and Computational Chemistry, Xiamen 361005, China*

---

## Repository Structure & Complete File Directory

```text
co_ru_2026/
├── Analysis/                 # Data post-processing and analysis toolkits
│   ├── analysis.txt          # Documentation for data processing workflows
│   ├── metadata3d.dat        # Structure logs for 3D coordinate sampling
│   ├── rmsdYxstat.cpp        # C++ routine to compute substrate Ru atom RMSD
│   ├── verticalangleud.cpp   # C++ routine for polar orientational distribution
│   ├── wham3dpara.c          # C-based 3D Weighted Histogram Analysis Method
│   └── wham3dpara.cpp        # C++ implementation of 3D-WHAM for PMF reconstruction
├── Dataset/                  # Raw training configurations for ML potentials
│   ├── 4x4x5_rucoint/        # Training configurations for CO-Ru(0001) chemisorption interaction
│   ├── 4x4x5_surf/           # Training geometries for the bare Ru(0001) substrate surface
│   ├── 4x4x5_threebody/      # Training geometries for CO-CO-Ru(0001) three-body lateral interaction
│   └── dpa2-finetune/        # Exhaustive data configs used for DPA-2 model fine-tuning
├── Examples/                 # Production-ready MD input configurations
│   ├── input-20205/          # US inputs for 0.500 ML coverage on a 20x20x5 Ru(0001) surface
│   ├── input-445/            # US inputs for 0.500 ML coverage on a 4x4x5 Ru(0001) surface
│   ├── metadata.dat          # Sampling spacing parameters for standard windows
│   └── metadata3d.dat        # Window indices for multi-dimensional enhanced sampling
├── MD_Code/                  # Master folder of molecular dynamics wrappers
│   ├── Umbrellasampling/     # Umbrella sampling wrappers for EANN & DPA-2 potentials
│   └── UnbiasedMD/           # Standard unbiased MD simulation codebase for CO diffusion
├── ML_Models/                # Machine-learning potential energy surface (ML-PES) interfaces
│   ├── 20x20x5_DPA2/         # Contains the deepmd frozen potential file `dpa2_m_OC20M247f.pth`
│   ├── 2x2x5_EANN/           # 6D rigid-surface potential and 42D full-flexible lattice EANN potentials
│   └── 4x4x5_EANN/           # High-dimensional EANN potential via Many-Body Expansion (MBE) strategy
├── UmbrellaSampling/         # Window-by-window implementation of Umbrella Sampling
│   ├── US-225/               # Code for 2x2x5 cell (INPUT for 6D PES, INPUT.2layer for 42D PES)
│   ├── US-445/               # Code for 4x4x5 cell (INPUT for two COs on a two-layer relaxed surface)
│   └── US-DPA2/              # Python driver script utilizing DeePMD-kit v3.1.1
└── UnbiasedMD/               # Unconstrained molecular dynamics trajectories
    ├── 2x2x5_EANN/           # Code for CO diffusion on a 2x2x5 Ru(0001) surface via EANN potential
    └── 4x4x5_EANN/           # Code for CO diffusion on a 4x4x5 Ru(0001) surface via EANN potential
```


## Detailed Directory Descriptions & Compilation Guide

### 1. Advanced Structural Analysis (`/Analysis`)

This section contains custom data post-processing tools written in high-performance C/C++ to analyze MD trajectories:

* **`wham3dpara.cpp` / `wham3dpara.c**`: Implements the multi-dimensional Weighted Histogram Analysis Method (WHAM) to eliminate biasing potentials and reconstruct full Potential of Mean Force (PMF) profiles.
* **`verticalangleud.cpp`**: Extracts the polar orientational angle $\alpha$ (the angle between the C–O bond vector and the surface normal) and applies a **Jacobian correction** (by dividing the raw frequencies by $\sin(\alpha)$) to isolate true orientational confinement $P(\alpha \mid z)$ from spherical geometric bias.
* **`rmsdYxstat.cpp`**: Quantifies the localized substrate relaxation and surface response induced by CO. It extracts the Root-Mean-Square Deviation (RMSD) of the top active Ru layers against an equilibrated bare Ru(0001) surface reference using a strict bonding cutoff of 2.20 Å.

### 2. Machine-Learning Potential Architectures (`/ML_Models`)

Trained neural-network potential parameter files and underlying evaluation wrappers:

* **`2x2x5_EANN/`**: Features the 6D analytical mapping potential (`nnfit_3.14.txt`, evaluated by `pes_ruFIXco.f90`) alongside the flexible 42D high-dimensional PES framework (`pes8.dat`, evaluated by `pes_hd.f90` / `pes_hd.inc1.f90`).
* **`4x4x5_EANN/`**: A complex potential constructed using the *Many-Body Expansion (MBE)* method to track lateral crowd interaction under high-coverage limits. It integrates cluster mapping (`pes_ococ_v4.2.F90`, `pes_ococ_w87.txt`), baseline matrices (`pes_pbe_Ru.f90`, `pes_tbody_int.f90`), and the comprehensive library `4x4x5_pes_coru`.
* **`20x20x5_DPA2/`**: Contains the model graph `dpa2_m_OC20M247f.pth` generated via active learning cycles to simulate large-scale catalysis pathways.

### 3. Molecular Dynamics Simulation Engines (`/UmbrellaSampling` & `/UnbiasedMD`)

Depending on the specific system setup, calculations are steered through separate compiler layers:

* **EANN Potentials (Fortran Core Environments)**:
Routines under `US-225`, `US-445`, and `UnbiasedMD` directories are automated through Makefiles. Move into the corresponding subfolder and compile directly:
```bash
make clean && make
./md_run.exe < INPUT

```


* For `US-445`, detailed execution profiles across all individual sampling steps can be synchronized using setup blocks stored inside `Examples/input-445/`.


* **DPA-2 Fine-Tuned Potentials (Python and DeePMD Core Frameworks)**:
Routines under `US-DPA2` handle ultra-large calculations via a pythonic pipeline:
```bash
cd UmbrellaSampling/US-DPA2
python dirofmain/main.py

```


* **Environment Requirement**: Execution requires a valid installation of **DeePMD-kit v3.1.1** configured with **CUDA 12.9** acceleration capability or higher.



---

## Training Data Availability (`/Dataset`)

To facilitate full methodology reproduction and deep-learning developments in heterogeneous on-surface catalysis, we share the exhaustive datasets used to fit the energy and force criteria:

* **`4x4x5_surf/`**: Snapshots tracking thermal vibrations and lattice flexibility of the clean Ru(0001) substrate.
* **`4x4x5_rucoint/`**: Local variations recording chemisorption configuration changes across atop, bridge, and threefold hollow chemical centers.
* **`4x4x5_threebody/`**: Structural data sets mapping out co-adsorbed repulsions and multi-molecule three-body terms.
* **`dpa2-finetune/`**: Global trajectory datasets compiled during active-learning iterations to fine-tune the macro-cell DPA-2 network.

```

```
