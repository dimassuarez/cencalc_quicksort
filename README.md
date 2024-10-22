# cencalc_quicksort

## Conformational ENtropy CALCulations using QuickSort

### Description

This version of CENCALC has been built by modifying the original code(https://github.com/ernestosuarez/CENCALC) developed by E. Suárez. The software has been designed to estimate the conformational entropy of single molecules from extended computer simulations, especially Molecular Dynamics (MD). On input CENCALC needs both trajectory coordinates and topology information in order to characterize the conformational states of the molecule of interest. The molecular conformers are identified by discretizing the time evolution of internal rotations. After this transformation, CENCALC determines the probability mass functions of the individual torsions and uses them for conformational entropy estimations. 

The version of CENCALC described here can use the classical Mutual Information Expansion (MIE) up to second-order, but the default method corresponds to the correlation-corrected multibody local approximation (CC-MLA). With respect to the original CENCALC code, the present version enables a faster and more exhaustive application of the CC-MLA technique, which can be combined with a distance-based cutoff criterion. In this case, CENCALC requires, as additional input, a distance matrix, for example, an interatomic distance matrix containing mean distance values derived from an MD trajectory in order to include only correlation effects among torsion angles whose mean separation is below a predefined cutoff. To speed up the CC-MLA calculations, an N x M matrix describing the evolution of M discretized dihedral angles along N MD snapshots is first transformed into a one-dimensional array of N big integer numbers. This array is efficiently sorted by using a parallel recursive implementation of the QUICKSORT algorithm, allowing thus to identify the different conformational states populated by the M dihedrals along the MD trajectory and to determine their relative abundances. Other methodologies for approaching to the full conformational entropy: the classical MIE at high orders, the Approximate MIE and the Multibody Local Approximation (MLA), are not available in this version and the original code should be used instead. 

The various techniques available in CENCALC have been discussed in the literature. Users are therefore encouraged to consult the references cited below before using the software. Note also that this manual is a revised version of that of the original code, which contains further practical information, and that the PDF included in the distribution gives also further details.   

We kindly request that any use of the CENCALC software or derivative should include at least the following citation:

1)	E. Suárez, N. Díaz, J. Méndez and D. Suárez. CENCALC: A Computational Tool for Conformational Entropy Calculations from Molecular Simulations. J. Comput. Chem. 54, 2031. DOI: 10.1002/jcc.23350

The methods implemented in CENCALC are fully described in the following references: 

2)	E. Suárez, N. Díaz and D. Suárez. Entropy Calculations of Single Molecules by Combining the Rigid-Rotor and Harmonic-Oscillator Approximations with Conformational Entropy Estimations from Molecular Dynamics Simulations J. Chem. Theor. Comput. 2011, 7, 2638-2653. DOI: 10.1021/ct200216n

3)	E. Suárez, D. Suárez. Multibody Local Approximation: Application to Conformational Entropy Calculations on Biomolecules. J. Chem. Phys. 2012, 137, 084115. DOI: 10.1063/1.4748104

4)	D. Suárez and N. Díaz. Toward Reliable and Insightful Entropy Calculations on Flexible Molecules. J. Chem. Theory Comput. 2022, 18, 7166–7178. DOI: 10.1021/acs.jctc.2c00858

All questions regarding the usage of this version of the CENCALC program or bug reports should be addressed to Dimas Suárez (dimas@uniovi.es). 

### Installation and Usage

The CENCALC software consists mainly of two standalone codes written in FORTRAN90, cencalc_prep.f90 and cencalc_ccmla.f90. The first program carries out various preparatory tasks prior to the main entropy calculations that are performed by cencalc_ccmla.f90. In this version, the cencalc_ccmla.f90 program makes use of the FORTRAN90 modules qsort.f90 , implementing the recursive quicksort algorithm as programmed by David Bal (https://bitbucket.org/daviddbal/multi-threaded-quicksort-mergesort-fortran/src/master/), and qsort_serial.f90, which implements the serial version. The FORTRAN90 codes take advantage of shared-memory parallel computers through the OpenMP Application Program Interface. Optionally, cencalc_prep.f90 can be linked to the  DISLIN high-level plotting library (https://www.dislin.de/) to make plots of probability mass functions and time evolution  of discretized dihedrals.  

One simple Makefile is provided to build the binaries using the GCC compilers (v. 8.3.1 or higher).

To facilitate the entropy estimations from MD trajectories generated by the AmberTools and/or Amber packages (http://ambermd.org), an auxiliary bash script (`calc_sconform.sh`) has been included. This script uses the cencalc_ccmla and cencalc_prep binaries as well as a separate python script (`get_tor.py`) that reads the topology information from AMBER parm files and identifies all the rotatable bonds that are required to characterize the conformational state of the molecule of interest. `calc_sconform.sh` may apply a composite scheme to calculate the total conformational entropy. It uses also a template python script (`sconform_plot.py`) to make plots of entropy curves. Other programs that are used by calc_sconform.sh are ncdump (from NetCDF) and GNU parallel. In principle, users could easily adapt `calc_sconform.sh` to work with other simulation packages and/or setup other specifically-tailored entropy calculations. 

### Example: Conformational Entropy of Amber Trajectories 

The example directory contains the required input topology and trajectory files to make a test run of the calc_sconform.sh script, which drives all the necessary steps and data transformations described in the previous section. The topology and trajectory files are located in the top_traj_files subdirectory. They correspond to a gas-phase 2.0 microsecond-long MD simulation of tenofovir using the GAFF force field. One million frames were saved for analysis. See the output_files subdirectory for cencalc output files.

Prior to its execution  , users must edit the `calc_sconform.sh` file  o set the variables pointing to the various programs being used by the script. These lines are placed at the top of the file and include the following variable definitions:

```
#CENCALC programs
export LD_LIBRARY_PATH="/opt/dislin/":$LD_LIBRARY_PATH # Optional
CENCALC_PATH="/home/dimas/SCRIPTS/CCMLA_QSORT/"

# AMBER trajectory analysis software
CPPTRAJ="$AMBERHOME/bin/cpptraj"

#  GNU Parallel builds and runs commands in parallel. 
PARALLEL="/opt/parallel-bash/bin/parallel --no-notice "  ```
```

Runtime options are passed to `calc_sconform.sh` as external variable declarations. For the processing of a MD trajectory, these options include the full topology filename, the directory containing the trajectory files, the alias (basename) of the trajectory files (with .mdcrd extension), the scratch directory, etc. 

For example, if `calc_sconform.sh` is copied into the example directory, the following command

```
env NPROCS=6 MOL_MASK=":1" REFTOP="$PWD/top_traj_files/drug_25.top" \
    TRAJDIR="$PWD/top_traj_files" PREFIX_MDCRD="drug_25" SUFFIX_MDCRD=".mdcrd" \
    DO_CC_MLA=1  CUTOFF=" -1 "  \
    SCRATCH="/scratch" ./calc_sconform.sh
```

will process the trajectory drug_25*.mdcrd file(s) located in the top_traj_files folder in order to obtain the time evolution of the discretized torsions and calculate the CC-MLA entropy (DO_CC_MLA=1) with no cutoff (CUTOFF=-1).  A modified `sconform_plot.py` script is also written in the output directory that produces  convergence plots of the entropy calculations. Successful execution of this test should create the following output files:

```
atoms_in_tor.info  d0002.dat.png  d0007.dat.png  d0012.dat.png  d0017.dat.png        MINIMA.dat               s_ccmla_-1.tab
ccmla_plot.tab     d0003.dat.png  d0008.dat.png  d0013.dat.png  d0018.dat.png        reduced_dist_matrix.dat  s_ccmla_plot.png
cencalc_prep.out   d0004.dat.png  d0009.dat.png  d0014.dat.png  distance_matrix.dat  s1.out                   s_ccmla.tab
CPU_TIME.dat       d0005.dat.png  d0010.dat.png  d0015.dat.png  matrix_cpptraj.out   s1_plot.tab              sconform_plot.py
d0001.dat.png      d0006.dat.png  d0011.dat.png  d0016.dat.png  MATRIX.dat           s_ccmla_-1.out           torsion.in
```
