# cencalc_quicksort

# Conformational ENtropy CALCulations using QuickSort

Description
This version of CENCALC has been built by modifying the original code developed by E. Suárez. The software has been designed to estimate the conformational entropy of single molecules from extended computer simulations, especially Molecular Dynamics (MD) simulations. On input CENCALC needs both trajectory coordinates and topology information in order to characterize the conformational states of the molecule of interest. The molecular conformers are identified by discretizing the time evolution of internal rotations. After this transformation, CENCALC determines the probability mass functions of the individual torsions and uses them for conformational entropy estimations. 

Although the version of CENCALC described here can use the classical Mutual Information Expansion (MIE) up to second-order, the default method corresponds to the correlation-corrected multibody local approximation (CC-MLA). With respect to the original CENCALC code, the present version enables a faster and more exhaustive application of the CC-MLA technique, which can be combined with a distance-based cutoff criterion. In this case, CENCALC requires, as additional input, a distance matrix, for example, an interatomic distance matrix containing mean distance values derived from an MD trajectory in order to include only correlation effects among torsion angles whose mean separation is below a predefined cutoff. To speed up the CC-MLA calculations, an N x M matrix describing the evolution of M discretized dihedral angles along N MD snapshots is first transformed into a one-dimensional array of N big integer numbers. This array is efficiently sorted by using a parallel recursive implementation of the QUICKSORT algorithm, allowing thus to identify the different conformational states populated by the M dihedrals along the MD trajectory and to determine their relative abundances. Other methodologies for approaching to the full conformational entropy: the classical Mutual Information Expansion (MIE) at high orders, the Approximate MIE (AMIE) and the Multibody Local Approximation (MLA), are not available in this version and the original code should be used instead. 

All the assumptions and equations defining the various techniques available in CENCALC have been discussed in the literature. Users are therefore encouraged to consult those references cited below before using the software. Note also that this manual is a revised version of that of the original code, which contains further practical information. 

We kindly request that any use of the CENCALC software or derivative should include at least the following citation:

1)	E. Suárez, N. Díaz, J. Méndez and D. Suárez. CENCALC: A Computational Tool for Conformational Entropy Calculations from Molecular Simulations. J. Comput. Chem. 54, 2031. DOI: 10.1002/jcc.23350.

The methods implemented in CENCALC are fully described in the following references: 

2)	E. Suárez, N. Díaz and D. Suárez. Entropy Calculations of Single Molecules by Combining the RigidRotor and Harmonic-Oscillator Approximations with Conformational Entropy Estimations from Molecular Dynamics Simulations J. Chem. Theor. Comput. 2011, 7, 2638-2653. DOI: 10.1021/ct200216n

3)	E. Suárez, D. Suárez. Multibody Local Approximation: Application to Conformational Entropy Calculations on Biomolecules. J. Chem. Phys. 2012, 137, 084115. DOI: 10.1063/1.4748104.

All questions regarding the usage of this version of the CENCALC program or bug reports should be addressed to Dimas Suárez (dimas@uniovi.es). 
