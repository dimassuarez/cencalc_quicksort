! This is modified version of the CENCALC_PREP code.
! The original version can be found at.
!
!  https://github.com/ernestosuarez/CENCALC
!
!===================================================================
!   CONFORMATIONAL ENTROPY CALCULATION
!
! Copyright (C) 2011 Ernesto Suarez Alvarez
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   http://www.gnu.org/licenses/>.
!-----------------------------------------------------------------------------------------
!  Changes of the original code written by Dimas Suarez (dimas@uniovi.es)  2022
!-----------------------------------------------------------------------------------------
!  Any use of the CENCALC software or derivative should include at least the following
!  citation:
!
!  1)E. Suarez, N. Diaz, J. Mendez and D. Suarez. CENCALC: A Computational Tool for
!    Conformational Entropy Calculations from Molecular Dynamics Simulations.
!    J Comput Chem.  2013 ;34(23):2041-54. doi: 10.1002/jcc.23350.
!
!  The methods implemented in CENCALC are fully described in the following references:
!
!  2)E. Suarez, N. Diaz and D. Suarez. Entropy Calculations of Single Molecules by
!    Combining the RigidRotor and Harmonic-Oscillator Approximations with Conformational
!    Entropy Estimations from Molecular Dynamics Simulations
!    J. Chem. Theor. Comput. 2011 , 7, 8, 2638-2653 doi: 10.1021/ct200216n
!
!  3)E. Suarez, D. Suarez. Multibody Local Approximation: Application in the Conformational
!    Entropy Calculation on Biomolecules.
!    J. Chem. Phys. 137, 084115 (2012);  doi: 10.1063/1.4748104
!
!-----------------------------------------------------------------------------------------

!     _____________________________
!    /|                           |
!    /|       MAIN PROGRAM        |
!    /|___________________________|
!    //////////////////////////////


!*****************************************************************************************
PROGRAM CENCALC_PREP
!*****************************************************************************************
! cencalc_prep
!
! This program prepares "MATRIX.dat" and "reduced_dist_matrix.dat"
! that are required to run centro_omp program.
!
! Files needed:
!
!        1-All the time series of each dihedral angle (one file by dihedral,
!            for example: file1.dat file2.dat ...). In each "file*.dat" the first
!            column is normally the time in (ps) that will not be read, and the
!            second which is the important one, is the dihedral angle value.
!            It also possible to read another column using the option-usecol.
!            For more info use the option-h/-help.
!
!
!        The Cut-off option is activated by default, in this case it
!        is also needed:
!
!        2-The distance matrix of the whole protein or system in vacuo
!            using the file distance_matrix.dat
!        3-The information about which atoms are involved in each dihedral.
!            The file atoms_in_tor.info should specify the two central atoms
!            of each angle. The first line corresponds to the first dihedral whose
!            time series is in file1.dat and so on.
!
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
! QUICK HELP:
!
! A quick help can be viewed from the command line by doing: centro_prep-help
!-----------------------------------------------------------------------------------------

   use parameters
   use omp_lib
   IMPLICIT NONE
   !-VARIABLE DEFINITIONS--------------------------------------------------------------
   integer NumFiles                                        !Number of data files
   integer UseCol                                          !Which column use in the data files
   integer NumAtoms                                        !Number of atoms in the distance matrix (lines or columns)
   integer NumSnap                                         !Number of snapshots
   integer NumProcs                                      !Number of Processors
   integer NumThreads                                    !Number of Threads
   character*60  DistMatrixName                            !Distance matrix file name
   character*60  InfoFileName                              !Info file name
   character*60, dimension (:), allocatable:: FileList      !List of data files (time evolution of torsions)
   character*60, dimension (:), allocatable:: IdDih         !Identification of Dih variables
   integer NumArg                                          !Number of arguments in the command line
   character*60, dimension (:), allocatable:: Arguments     !Arguments in the command line
   character*3 simplify                                    !Eliminate conformationally rigid torsions? (yes/no)
   logical UseCut                                          !Use cutoff? (True/False)
   logical RdBigMat                                        ! Read Big Matrix and Jump to Statistical Part
   logical WrBigMat                                        ! Write Big Matrixa and Stops
   character*120  BigMatFile                                ! BigMatirx filename
   integer, dimension (:,:), allocatable:: BigMatrix  ! Discretized data matrix before removing rigid torsions
   integer, dimension (:,:), allocatable:: Matrix    !Discretized data matrix after removing rigid torsions
   real(DP), dimension(:,:), allocatable  :: Density         !Von Mises PDF for plotting purposes
   real(DP), dimension(:,:), allocatable  :: HistDensity     !Histogram PDF for plotting purposes
   integer  plotsize                                       ! Integer scaling actor for plot size ( 1, 2, 3 ...)
   integer, dimension (:,:), allocatable:: Atoms_in_Dih     !Two-column matrix where is saved the content of InfoFileName
   real(DP), dimension (:,:), allocatable:: Dist            !Distance matrix
   real(DP), dimension (:,:), allocatable:: NewDistMat      !Output: reduced distance matrix
   integer NewNumCol                                       !New number of columns (the columns in reduced dist. matrix)
   real(DP), dimension (:,:), allocatable:: ang             !Dihedral angle value in degrees
   integer, dimension (:,:), allocatable:: DiscreteAng  ! Discretized dihedral angle
   real(DP) k_value                                        !(see explanation in default values)
   logical  AnalyticGrad                                   !Use analytic gradient in the optimization?
   logical  Plot, Saveplotdata                            !Plot PDF data and time evolution for each dihedral ?
   real(DP) Step                                           !Step for the steepest-descendent optimization
   integer  MaxIterations                                  !Max number of iteration in the optimization
   real(DP) ConvCriterion                                  !Convergence criterion
   integer  MaxNumConf                                     !Max Number of Conformers to look for
   integer, dimension(:), allocatable ::   FNumMins        !final num of minimums
   real(DP), dimension(:,:), allocatable:: minimum         !Positions of the minimums in the PDF
   real(DP), dimension(:,:), allocatable:: Pop             !Population of each conformational state
   real(DP), dimension(:,:), allocatable:: tlife           !Average life (in # snap) of discretized conformational states
   integer, dimension(:,:), allocatable:: nlife            !Number of consecutive appearances of conformational states
   integer, dimension(:,:), allocatable:: maxlife          !Maximum life (in # snap) of discretized conformational states
   integer, dimension(:,:), allocatable:: minlife          !Minimum life (in # snap) of discretized conformational states (Currently, without use)
   integer, dimension(:), allocatable:: ilife              ! ilife(istate) line in # of consecutive snap of a conformational state
   integer, dimension(:), allocatable:: itorlife           ! Dihed ID associated to ilife(istate)
   real(DP), parameter:: LimTOR = 0.0005d0                     ! Lower limit for the ratio ilife(istate)/Maxlife for a torsional state to be considered as such for statistical analysis
   integer, dimension(:), allocatable:: ifreq              ! Number of discretized conformational changes/Pointer
   real(DP), dimension(:), allocatable:: freq              ! Frequency of discretized conformational changes
   real(DP), dimension(:), allocatable:: rtmpaux          ! Frequency of discretized conformational changes
   real(DP), dimension(:), allocatable:: logfreq           ! Log10-Frequency of discretized conformational changes
   real(DP), dimension(:), allocatable:: s1                !First order entropy of dihedral angles
   real(DP) stot,sfree                                           !Total First order entropy
   real(DP), parameter:: s1_limit = 0.01d0                     ! Minimum entropy value for considering a dihedral as conformationally active.
   real(DP), parameter:: freq_fast = 0.10d0                    ! Minimum rate of conformational change per #snap of "fast" dihedrals
   real(DP), parameter:: freq_medium = 0.001d0                 ! Minimum rate of conformational change per #snap of "not too fast, but not too slow" dihedrals
   integer  nfast, nmedium, nslow, nquasifrozen
   integer  nmin
   integer, parameter:: maxgroup= 20                    ! Minimum rate of conformational change per #snap of "not too fast, but not too slow" dihedrals
   integer, dimension(:), allocatable:: quasifrozen        ! pointer to quasi frozen dihedrals
   integer, dimension(:), allocatable:: fast               ! pointer to fast dihedrals
   integer, dimension(:), allocatable:: medium             ! pointer to medium dihedrals
   integer, dimension(:), allocatable:: slow               ! pointer to slow dihedrals
   integer, dimension(:), allocatable:: active             ! pointer to selected active dihedrals
   integer, dimension(:,:), allocatable  :: corr_group     ! Von Mises PDF for plotting purposes
   integer, dimension(:), allocatable:: ngroup             ! auxiliary array
   real CpuTimeTot                                         !CPU time variables
   integer(K4C) RTime1,RTime2,rate                 !Auxiliary variables to compute real-time
   integer    lgroup(maxgroup)
   character*12, dimension(:), allocatable:: label_freq    ! Text Label of dihedrals
   real(DP) dummy, tsum ,rinf, rsup, freq_thres, w         !
   character*1 charact                                     ! Auxiliary or Dummy variables
   character*120 cdummy                                     !
   logical minima_file
   integer idummy, idummy2, iargc, ios , ig, ng, logfreq_min, logfreq_max                   !
   integer i, j, k, l, ithread, C, A, iang                        !
   integer ii, jj, kk, ll, istate, nstate                       !
   !-----------------------------------------------------------------------------------

!  

   write(*,'(A33)') " _______________________________ "
   write(*,'(A33)') "||  _____                       |"
   write(*,'(A33)') "||  \\   ||                     |"
   write(*,'(A33)') "||   \\  PROGRAM CENCALC_PREP   |"
   write(*,'(A33)') "||   //          v0.3           |"
   write(*,'(A33)') "||  //___||                     |"
   write(*,'(A33)') "||______________________________|"
   print*," "

   call system_clock(count=RTime1, count_rate=rate)
!
!---DEFAULT VALUES--------------------------
   UseCol = 2                                                !Which column is going to be discretized?

   simplify="yes"                                          !Eliminate columns with null entropy

   InfoFileName="atoms_in_tor.info"                        !Which atoms define de dihedral

   DistMatrixName="distance_matrix.dat"                    !All solute atoms distance matrix

   UseCut=.true.                                           !Use cut-off

   k_value = 0.50D0                                           !The smoothing parameter "v" of
   !the von-Mises kernel density estimation
   !depends of the k_value see eq.(7) in ref:
   !Computational Statistics & Data Analysis
   !Volume 52, Issue 7, 15 March 2008, Pages 3493-3500.
   !Here we do not estimate k_value, but simply
   !set the value empirically in order to
   !slightly over-smooth the PDFs. In any case, the user
   !can change this value by the option "-k".

   MaxIterations = 1000                                      !Max number of Iterations in order to
   !find the PDFs minimums
   Step = 5.d0                                                  !Step-size during the optimization

   ConvCriterion = 1D-4                                      !Convergence criterion (gradient)

   MaxNumConf = 4                                          !Max number of conformers by torsion to look for

   AnalyticGrad=.false.                                    !Use analytical gradient for optimization,
   !the ".false." means the usage of linear
   !gradient interpolation
   Plot=.false.                                            !Plot Data
   Saveplotdata=.false.                                    ! Save Plot Data (may be useful for getting better plots using other software)

   RdBigMat=.false.
   WrBigMat=.false.
   BigMatFile='BIGMATRIX.dat'
   minima_file=.false.
!------------------------------------------

   NumArg = iargc()
   allocate(Arguments(NumArg))
   allocate(FileList(NumArg))
   allocate(IdDih(NumArg))

!---READING OPTIONS----------
   CALL Read_Options(FileList, UseCol, InfoFileName, DistMatrixName, &
      NumArg, NumFiles, simplify, UseCut, k_value, AnalyticGrad, Plot, Saveplotdata, Step, &
      ConvCriterion, MaxIterations, MaxNumConf, RdBigMat, WrBigMat, BigMatFile)


   print*," "
   print*," OPTIONS:"
   print*,"--------------------------------------------------"
   if(UseCut) then
      write(*,'(A32)') "Distance Matrix file name:     "
      write(*,*) DistMatrixName
      write(*,'(A32)') "Dihedral Information file name:"
      write(*,*) InfoFileName
      write(*,'(A25, A7)') "Using CutOff:           ","YES"
   else
      write(*,'(A25, A7)') "Using CutOff:           ","NO"
   endif
   write(*,'(A25, I7)') "Using column:           ",UseCol
   write(*,'(A25, F7.2)') "k_value:                ",k_value
   write(*,'(A25, I7)') "Max. num. of Iterations:",MaxIterations
   write(*,'(A25, F7.2)') "Step:                   ",Step
   write(*,'(A25, F7.2)') "Convergence Criterion:  ",ConvCriterion
   write(*,'(A25, I7)') "Max. num. of Conformers:",MaxNumConf
   if(AnalyticGrad)      write(*,'(A25, A7)') "Using Analytic Gradient:","YES"
   if(.not.AnalyticGrad) write(*,'(A25, A7)') "Using Analytic Gradient:","NO"
   write(*,'(A32)') "NumFiles:   "
   write(*,*) NumFiles
   print*," "
   print*, " FILES:"

   A = 0
   do i = 1, NumFiles
      k = i+6
      open(unit = k, file = FileList(i), status="OLD",&        !Just to check if the file exists
         iostat = ios)
      if(ios .ne. 0) then
         print*," "
         write(*,'(2A22)') "ERROR during opening: ",FileList(i)
         stop
      endif
      ios = 0
      IdDih(i)='NULL'
      read(k, '(A)') cdummy
      if ( cdummy(1:7) == '#Frame ')  then
         l = len(cdummy)
         do while ((cdummy(l:l).ne.' ') .and. (l .gt. 0))
            l = l-1
         enddo
         IdDih(i)=cdummy(8:l)
      else
         backspace(k)
      endif
      if(i .eq. 1) then
         do while(ios == 0)
            read(k, '(A1)',iostat = ios) charact          !Just to know the number of snapshots
            if(ios > 0) then
               print*, "Error(0) while reading ",FileList(i)
               stop
            endif
            if(ios == -1) exit
            A = A+1
         enddo
      endif
      if ( IdDiH(i) .ne. 'NULL') then
         write(*,'(I5, 1X, A10, A60)') i, FileList(i), IdDih(i)
      else
         write(*,'(I5, 1X, A10)') i, FileList(i)
      endif
      close(k)
   enddo
   NumSnap = A
   print*,"--------------------------------------------------"
   print*," "
   print*,"Number of Snapshots: ",NumSnap
   print*,"Info: The Number of Snapshots is read from the first file: ",FileList(1)
   print*," "
!---------------------------------------------------------------------------------
   if(UseCut) then
      open(unit = 9, file = DistMatrixName, status="OLD",iostat = ios)
      if(ios .ne. 0) then
         print*," "
         write(*,'(2A22)') "ERROR during opening: ",DistMatrixName
         stop
      endif
      close(9)
   endif
!---------------------------
   NumProcs = 1
   NumThreads = 1
!$omp parallel
!$ NumProcs = OMP_get_num_procs()
!$ NumThreads = OMP_get_num_threads()
!$omp end parallel

! Number of Threads cannot be greater than NumProcs since
! we request paralell IO of data files
!$ print*,'OMP Parallel Execution'
!$ print*,"Number of Processors=",NumProcs
!$ print*,"Number of Threads   =",NumThreads
!$ if ( NumThreads > NumProcs ) then
!$    call OMP_set_num_threads(NumProcs)
!$    print*,'Number of Threads set to ',NumProcs
!$ endif


!  Allocating Arrays
   allocate(BigMatrix(NumSnap, NumFiles))
   allocate(FNumMins(NumFiles))
   allocate(minimum(9, NumFiles))
   allocate(Density(0:360, NumFiles))
   allocate(Pop(9, NumFiles))
   allocate(tlife(9, NumFiles))
   allocate(nlife(9, NumFiles))
   allocate(maxlife(9, NumFiles))
   allocate(minlife(9, NumFiles))
   allocate(ilife(NumSnap))
   allocate(itorlife(NumSnap))
   allocate(ifreq(NumFiles))
   allocate(quasifrozen(NumFiles))
   allocate(fast(NumFiles))
   allocate(medium(NumFiles))
   allocate(slow(NumFiles))
   allocate(active(NumFiles))
   allocate(freq(NumFiles))
   allocate(rtmpaux(NumFiles))
   allocate(logfreq(NumFiles))
   allocate(ngroup(maxgroup))
   allocate(corr_group(maxgroup,NumFiles))
   allocate(label_freq(NumFiles))
   allocate(s1(NumFiles))
   allocate(HistDensity(0:72, NumFiles))
   allocate(ang(NumSnap, NumProcs))
   allocate(DiscreteAng(NumSnap, NumProcs))
!
!
!
   if ( RdBigMat ) then
      write(*,'(A25,1X,A60)') "Reading BIGMATRIX in ",BigMatFile
      open(12, file=BigMatFile,status="old")
      do i = 1, NumSnap
         read(12, '(10000I2)') (Bigmatrix(i, k), k = 1, NumFiles)
      enddo
      close(12)
      print*,"OK"
! Computing Pop Matrix
      do k = 1, NumFiles
         do i = 1, 9
            Pop(i, k)=0.0d0
         enddo
      enddo
!$omp parallel shared(BigMatrix,NumSnap,Pop) private(i,j,k)
!$omp do schedule(static, 1)
      do k = 1, NumFiles
         do i = 1, NumSnap
            j = BigMatrix(i, k)
            Pop(j, k)=Pop(j, k)+1.0d0
         enddo
      enddo
!$omp end do
!$omp end parallel
      do k = 1, NumFiles
         do ii = 1, 9
            Pop(ii, k)=100.0d0*(Pop(ii, k)/dfloat(NumSnap))
         enddo
         FNumMins(k)=maxval(BigMatrix(:,k))-minval(BigMatrix(:,k))+1
      enddo
      inquire(file='MINIMA.dat',Exist=minima_file)
      if ( minima_file) then
         open (12,file='MINIMA.dat',status='old')
         do k=1,NumFiles
            read(12,'(34X,9F12.2)',end=1111,err=1111) (minimum(i, k), i = 1, FNumMins(k))
         enddo
         close(12)
      endif

      GOTO 1111
   endif

!--CREATING THE BIG MATRIX---------
   write(*,'(A49)') " Codifying and creating the full data matrix..."
   do i = 1, NumFiles, NumThreads
      if ( i+NumThreads-1 .gt. NumFiles) then
         j = NumFiles
      else
         j = i+NumThreads-1
      endif
!
!$omp parallel default(shared) private(ithread, A, ios, dummy, k, kk)
!$omp do schedule(dynamic, 1)
!
      do ithread = 1, j-i+1
         k = i+ithread-1
         open(unit = k+6, file = FileList(k), status="OLD")
         if ( IdDiH(k) .ne. 'NULL') read(k+6,*)
         A = 1
         ios = 0
         do while((ios == 0) .and. ( A .lt. NumSnap) )
            read(k+6, *,iostat = ios) (dummy, kk = 1, (UseCol-1)), ang(A, ithread)
            do while((ang(A, ithread).lt.0.d0).or.(ang(A, ithread).ge.360.d0))
               if(ang(A, ithread).lt.0.d0) ang(A, ithread)=ang(A, ithread)+360.d0
               if(ang(A, ithread).ge.360.d0) ang(A, ithread)=ang(A, ithread)-360.d0
            enddo
            if(ios > 0) then
               print*, "Error while reading ",FileList(k)
               stop
            endif
            if(ios == -1) exit
            A=A+1
         enddo
         if(NumSnap .ne. A) then
            print*," "
            print*, "ERROR: The files ",FileList(1), " and ",FileList(k), &
               " have different number of snapshots or lines"
            stop
         endif
         close(k+6)
      enddo
!
!$omp end do
!$omp end parallel
!

!$omp parallel default(shared) private(ithread, k)
!$omp do schedule(dynamic, 1)
      do ithread = 1, j-i+1
         k = i+ithread-1
         call Get_PDF_Min(ang(:,ithread), NumSnap, k_value,    &
            AnalyticGrad, Step, ConvCriterion, MaxIterations,   &
            MaxNumConf, minimum(:,k), FNumMins(k), &
            Density(0:360, k), HistDensity(0:72, k))
         write(*,'(A20, ''# Minima='',I3, 2X, 9F12.2)') FileList(k), &
            FNumMins(k), (minimum(ii, k), ii = 1, FNumMins(k))
      enddo
!$omp end do
!$omp end parallel


!$omp parallel default(shared) private(ithread, k, ii, jj)
!$omp do schedule(static, 1)
      do ithread = 1, j-i+1
         k = i+ithread-1
         if ( FNumMins(k).eq.1) then
            do ii = 1, NumSnap
               DiscreteAng(ii, ithread)=1
            enddo
         else
            do ii = 1, NumSnap
               do jj = 1, FNumMins(k)
                  if (jj .eq. 1) then
                     if((ang(ii, ithread).lt.minimum(jj, k)).or.    &
                        (ang(ii, ithread).ge.minimum(FNumMins(k), k))) DiscreteAng(ii, ithread)=jj
                  else
                     if((ang(ii, ithread).lt.minimum(jj, k)).and.    &
                        (ang(ii, ithread).ge.minimum(jj-1, k))) DiscreteAng(ii, ithread)=jj
                  endif
               enddo
            enddo
         endif
      enddo
!$omp end do
!$omp end parallel


!$omp parallel default(shared) private(ithread, k, ii)
!$omp do schedule(static, 1)
      do ithread = 1, j-i+1
         k = i+ithread-1
         do ii = 1, NumSnap
            BigMatrix(ii, k)=DiscreteAng(ii, ithread)
         enddo
      enddo
!$omp end do
!$omp end parallel
!
   enddo
   print*,"OK"
   print*,"  "
!----------------------------------
! COMPUTING POPULATION OF EACH CONF STATE IN BigMatrix
   do k = 1, NumFiles
      do ii = 1, 9
         Pop(ii, k)=0.0d0
      enddo
   enddo
   do i = 1, NumFiles, NumThreads
      if ( i+NumThreads-1 .gt. NumFiles) then
         j = NumFiles
      else
         j = i+NumThreads-1
      endif
!$omp parallel default(shared) private(ithread, k, ii, jj)
!$omp do schedule(static, 1)
      do ithread = 1, j-i+1
         k = i+ithread-1
         do ii = 1, NumSnap
            jj = BigMatrix(ii, k)
            Pop(jj, k)=Pop(jj, k)+1.0d0
         enddo
      enddo
!$omp end do
!$omp end parallel
   enddo
   do k = 1, NumFiles
      do ii = 1, 9
         Pop(ii, k)=100.0d0*(Pop(ii, k)/dfloat(NumSnap))
      enddo
   enddo
!
! PLOTTING info
!
   if ( Plot ) then
      plotsize = 2
      write(*,'(A49)') " Creating PDF and Time evolution plots  ..."
      do k = 1, NumFiles
         call HistPDF_Time_Plot(Saveplotdata, NumSnap, Filelist(k), IdDiH(k),  &
            FNumMins(k), minimum(:,k), Pop(:,k),                   &
            Bigmatrix(:,k), &
            Density(0:360, k), HistDensity(0:72, k), plotsize)
      enddo
      print*,"OK"
      print*,"  "
   endif
!
!     Writing BigMatrix if requested and then stops
!
   if ( WrBigMat   ) then
      write(*,'(A25,1X,A60)') " Writing BIGMATRIX in ",BigMatFile
      open(12, file=BigMatFile,status='unknown')
      do i = 1, NumSnap
         write(12, '(10000I2)') (Bigmatrix(i, j), j = 1, NumFiles)
      enddo
      close(12)
      call cpu_time(CpuTimeTot)
      call system_clock(count=RTime2)

      write(*,'(A26,F11.2,A8,F9.2,A6)')" CPU-TIME : ",&
      CpuTimeTot," seconds",(CpuTimeTot/3600.)," hours"
      write(*,'(A26,F11.2,A8,F9.2,A6)')" REAL_TIME : ",&
      real(RTime2-RTime1)/rate," seconds",(real(RTime2-RTime1)/(rate*3600.))," hours"

      stop
   endif
!
! COMPUTING STATISTICS OF DIHEDRAL DYNAMICS
!
1111 CONTINUE
!
   do k = 1, NumFiles
      ifreq(k)=0
      freq(k)=0.0d0
      rtmpaux(k)=0.0d0
   enddo
   do k = 1, NumFiles
      do ii = 1, 9
         tlife(ii, k)=0.0d0
         nlife(ii, k)=0
         maxlife(ii, k)=0
         minlife(ii, k)=999999
      enddo
   enddo
!
   do k = 1, NumFiles
      iang = BigMatrix(1, k)
      do ii = 1, NumSnap
         ilife(ii)=0
         itorlife(ii)=0
      enddo
      istate = 1
      do ii = 2, NumSnap-1
         if ( iang .ne. BigMatrix(ii, k) ) then
            ilife(istate)=ilife(istate)+1
            itorlife(istate)=iang
            iang = BigMatrix(ii, k)
            istate = istate+1
         else
            ilife(istate)=ilife(istate)+1
         endif
      enddo
      if ( iang .eq. BigMatrix(NumSnap, k) ) then
         ilife(istate)=ilife(istate)+1
         itorlife(istate)=iang
      else
         ilife(istate)=ilife(istate)+1
         itorlife(istate)=iang
         istate = istate+1
         ilife(istate)=ilife(istate)+1
         itorlife(istate)=BigMatrix(NumSnap, k)
      endif
      nstate = istate
      do istate = 1, nstate
         iang = itorlife(istate)
         if (maxlife(iang, k).lt.ilife(istate))  maxlife(iang, k)=ilife(istate)
         if (minlife(iang, k).gt.ilife(istate))  minlife(iang, k)=ilife(istate)
      enddo
      do istate = 1, nstate
         iang = itorlife(istate)
         if (dfloat(ilife(istate))/dfloat(maxlife(iang, k)).gt.LimTOR) then
            nlife(iang, k)=nlife(iang, k)+1
            tlife(iang, k)=tlife(iang, k)+dfloat(ilife(istate))
            ifreq(k)=ifreq(k)+1
         endif
      enddo
   enddo
!
   do k = 1, NumFiles
      freq(k)=dfloat(ifreq(k))/dfloat(NumSnap)
      rtmpaux(k)=freq(k)
      ifreq(k)=k
      label_freq(k)=' '
      do ii = 1, 9
         if (nlife(ii, k) .gt. 0) tlife(ii, k)=tlife(ii, k)/dfloat(nlife(ii, k))
      enddo
   enddo
!
!     First order entropy is readily computed and allows the rapid
!     detection of frozen and quasi frozen dihedrals
!
   stot = 0.0d0
   nquasifrozen = 0
   do k = 1, NumFiles
      s1(k)=0.0d0
      do ii = 1, FNumMins(k)
         if ((Pop(ii, k)/100.0d0)  .gt. 0.0001d0 )  &
            s1(k)=s1(k)+(Pop(ii, k)/100.0d0)*log(Pop(ii, k)/100.0d0)
      enddo
      s1(k)=-R*s1(k)
      if (( s1(k) .lt.  s1_limit ) .and. ( FNumMins(k) .le. 1)) then
         label_freq(k)=' FROZEN '
      else if ( s1(k) .lt.  s1_limit ) then
         label_freq(k)=' QUASIFROZEN'
         nquasifrozen = nquasifrozen+1
         quasifrozen(nquasifrozen)=k
      endif
      stot = stot+s1(k)
   enddo
!
!     Prior to printing out information dihedral angles are sorted in terms of their
!     rate of conformational change and then categorized

   if (Numfiles .gt. 1) call sort2(NumFiles, rtmpaux, ifreq)
!
!    The rest of the dihedrals are classified as fast-medium-slow
!    dihedrals referring to their rate of conf. change per #snap-1
!
   nfast = 0
   nmedium = 0
   nslow = 0
   do l = 1, NumFiles
      k = ifreq(l)
      if (index(label_freq(k), 'FROZEN') .eq. 0) then
         if ( freq(k) .gt. freq_fast ) then
            label_freq(k)=' FAST '
            nfast = nfast+1
            fast(nfast)=k
         else if ( freq(k) .gt. freq_medium ) then
            label_freq(k)=' MEDIUM '
            nmedium = nmedium+1
            medium(nmedium)=k
         else
            label_freq(k)=' SLOW '
            nslow = nslow+1
            slow(nslow)=k
         endif
      endif
   enddo
!
!
!
   write(6, '(''%======================================================='')')
   write(6, '(''% Statistical Analyses of Discretized Dihedral Dynamics'')')
   write(6, '(''%======================================================='')')
   write(6, '(''%'',19X, ''Dihed ID  '',30X, ''S1(cal/mol K)  PDF minima (deg) '',&
      3X, ''% Population '',8X,  &
      ''Average Life (#snap)'',5X, ''Freq of Change (1/#snap) '')')
   do l = 1, NumFiles
      k = ifreq(l)

      if ( FNumMins(k) .le. 1 ) then
         write(6, '(''% '',A9, A50, 1X, F9.3, 2X, 1(1X, F5.1), 2(2X, ''--'',2X), 1(1X, F5.1), &
            2(2X, ''--'',2X), 1(F10.1), 2(4X, ''--'',4X), F9.6, A12)') &
            Filelist(k), IdDiH(k), s1(k), (minimum(ii, k), ii = 1, FNumMins(k)), (Pop(ii, k), ii = 1, FnumMins(k)) &
            ,(tlife(ii, k), ii = 1, FNumMins(k)), freq(k), label_freq(k)
      elseif ( FNumMins(k) .eq. 2 ) then
         write(6, '(''% '',A9, A50, 1X, F9.3, 2X, 2(1X, F5.1), 2X, ''--'',2X, &
            2(1X, F5.1), 2X, ''-'',2X, 2(F10.1), 4X, ''--'',4X, F9.6, A12)') &
            Filelist(k), IdDiH(k), s1(k), (minimum(ii, k), ii = 1, FNumMins(k)), (Pop(ii, k), ii = 1, FnumMins(k)) &
            ,(tlife(ii, k), ii = 1, FNumMins(k)), freq(k), label_freq(k)
      elseif ( FNumMins(k) .eq. 3 ) then
         write(6, '(''% '',A9, A50, 1X, F9.3, 2X, 6(1X, F5.1), 3(F10.1), F9.6, A12)') &
            Filelist(k), IdDiH(k), s1(k), (minimum(ii, k), ii = 1, FNumMins(k)), (Pop(ii, k), ii = 1, FnumMins(k)) &
            ,(tlife(ii, k), ii = 1, FNumMins(k)), freq(k), label_freq(k)
      endif
   enddo
   write(6, '(''% Lower freq limit for FAST dihedrals     = '',F12.9)') freq_fast
   write(6, '(''% Lower freq limit for MEDIUM dihedrals   = '',F12.9)') freq_medium
   write(6, '(''% Maximum S1 for QUASIFROZEN dihedrals    = '',F12.3)') s1_limit
   write(6, '(''% Total first order entropy (cal/(mol K)) = '',F12.3)') stot
   write(6, '(''%======================================================='')')


!--SELECTING THOSE COLUMNS OF BigMatrix THAT CHANGE CONFORMATIONALLY

   write(6,'(A)') "% First selection of active columns (torsions)."
   C = 0
   if ( simplify .ne. 'yes' ) then
      do j = 1, nquasifrozen
         C = C+1
         active(C)=quasifrozen(j)
      enddo
      write(6,'(A)') "% Removing FROZEN dihedrals "
   else
      write(6,'(A)') "% Removing both FROZEN and QUASIFROZEN dihedrals "
   endif
   do j = 1, nslow
      C = C+1
      active(C)=slow(j)
   enddo
   do j = 1, nmedium
      C = C+1
      active(C)=medium(j)
   enddo
   do j = 1, nfast
      C = C+1
      active(C)=fast(j)
   enddo
   NewNumCol = C             !New number of columns after the elimination
   write(6,'(''% New number of columns (dihedrals) in reduced data matrix='',I4)') NewNumCol
   C = 1
   if ( ( simplify .ne. 'yes' ) .and. (nquasifrozen.gt.0)) then
      print*,"% quasifrozen dihedrals=",nquasifrozen, " Columns=",C, nquasifrozen
      C = C+nquasifrozen
   endif
   if ( nslow .gt. 0 ) then
      write(6,'(''% Slow dihedrals='',I3, '' Columns= '',I3,''  '',I3)') nslow, C, C+nslow-1
      C = C+nslow
   endif
   if ( nmedium .gt. 0 ) then
      write(6,'(''% Medium-rate dihedrals='',I3, '' Columns= '',I3,''  '',I3)') nmedium, C, C+nmedium-1
      C = C+nmedium
   endif
   if ( nfast .gt. 0 ) then
      write(6,'(''% Fast dihedrals='',I3, '' Columns= '',I3,''  '',I3)') nfast, C, C+nfast-1
      C = C+nfast
   endif
   print*,"  "

!
!--CREATING THE "SMALL" MATRIX
!
   allocate(Matrix(NumSnap, NewNumCol))
   do j = 1, NewNumCol
      do k = 1, NumSnap
         Matrix(k, j)=BigMatrix(k, active(j))
      enddo
   enddo
   print*,"  "
!--------------------------------------------------------------

!----OPENING AND READING THE DISTANCE MATRIX-----------
   if(UseCut) then
      write(*,'(A33)')" Reading the distance matrix ... "
      open(unit = 9, file = DistMatrixName, status="OLD",iostat = ios)
      NumAtoms = 0
      ios = 0
      do while(ios == 0)
         read(9, *,iostat = ios) dummy                         !Just to know the number of lines (ATOMS)
         if(ios > 0) then
            print*, "Error while reading ",DistMatrixName
            stop
         endif
         if(ios == -1) exit
         NumAtoms = NumAtoms+1
      enddo
      rewind(9)
      allocate(Dist(NumAtoms, NumAtoms))

      DO i = 1, NumAtoms
         READ(9, *)(Dist(i, j), j = 1, NumAtoms)
      ENDDO
      close(9)
      print*,"OK"
      print*,"  "
   endif
!-----------------------------------------------------


!----OPENING AN READING THE INFO FILE-----
!
   if(UseCut) then
      write(*,'(A27)')" Reading the info file ... "

      allocate(Atoms_in_Dih(NumFiles, 2))
      open(unit = 10, file = InfoFileName, status="OLD",iostat = ios)
      if(ios .ne. 0) then
         print*," "
         write(*,'(2A22)') "ERROR during opening: ",InfoFileName
         stop
      endif

      c = 0
      do while(ios == 0)
         read(10, *,iostat = ios) idummy                       !Just to know the number of lines
         if(ios > 0) then
            print*, "Error while reading ",InfoFileName
            stop
         endif
         if(ios == -1) exit
         c = c+1
      enddo
      rewind(10)
      if(c .gt. NumFiles) then
         print*, "WARNING: The number of lines in the info file is greater"
         print*, "         than the number of dihedrals (files)"
         print*, " "
      elseif(c .lt. NumFiles) then
         print*, "ERROR: The number of dihedrals (files) is greater than"
         print*, "       the number of lines in the info file, the program"
         print*, "       can not countinue using cut-off."
         print*, " "
         stop
      endif

      c = 0
!
! This must be recoded by reading all dihedral info!!!
!
      do i = 1, NumFiles
         read(10, *,iostat = ios) (Atoms_in_Dih(i, j), j = 1, 2)
         if(ios .ne. 0) then
            print*, "Error while reading ",InfoFileName
            stop
         endif
      enddo
!
      close(10)
      print*,"OK"
      print*,"  "
   endif
!---------------------------------------


!---MAKING THE NEW (OR REDUCED) DISTANCE MATRIX---
   if(UseCut) then
      write(*,'(A49)')" Making the new (or reduced) distance matrix ... "
      allocate(NewDistMat(NewNumCol, NewNumCol))
      do i = 1, NewNumCol
         do j = 1, NewNumCol
            ii = active(i)
            jj = active(j)
            if(i .eq. j) then
               NewDistMat(i, j)=0.0d0
            else
               NewDistMat(i, j)=(0.25d0)*( Dist(Atoms_in_Dih(ii, 1), Atoms_in_Dih(jj, 1)) &
                  + Dist(Atoms_in_Dih(ii, 2), Atoms_in_Dih(jj, 2)) &
                  + Dist(Atoms_in_Dih(ii, 1), Atoms_in_Dih(jj, 2)) &
                  + Dist(Atoms_in_Dih(ii, 2), Atoms_in_Dih(jj, 1)))
            endif
         enddo
      enddo
      print*,"OK"
      print*,"  "
   endif
!-------------------------------------------------


!----PRINTING OUTPUT------
   !writing the reduced distance matrix
   if(UseCut) then
      write(*,'(A50)')" Writing in reduced_dist_matrix.dat ... "
      open(11, file="reduced_dist_matrix.dat")

      do I = 1, NewNumCol
         write(11, '(10000F9.3)')(NewDistMat(i, j), j = 1, NewNumCol)
      enddo
      close(11)
      print*,"OK"
      print*,"  "
   endif

! Writing the final matrix
   write(*,'(A26)')" Writing in MATRIX.dat ... "
   open(12, file="MATRIX.dat")
   do i = 1, NumSnap
      write(12, '(10000I2)') (Matrix(i, j), j = 1, NewNumCol)
   enddo
   close(12)

   call cpu_time(CpuTimeTot)
   call system_clock(count=RTime2)

   write(*,'(A26,F11.2,A8,F9.2,A6)')" CPU-TIME : ",&
      CpuTimeTot," seconds",(CpuTimeTot/3600.)," hours"
   write(*,'(A26,F11.2,A8,F9.2,A6)')" REAL_TIME : ",&
      real(RTime2-RTime1)/rate," seconds",(real(RTime2-RTime1)/(rate*3600.))," hours"


   print*,"OK"
   print*,"  "
!--------------------------

   STOP " DONE!"
!----------------------------------------------------------------------------------------
END PROGRAM
!----------------------------------------------------------------------------------------

!      _____________________________
!     /|                           |
!     /| SUBROUTINES AND FUNCTIONS |
!     /|___________________________|
!     //////////////////////////////


!*****************************************************************************************
SUBROUTINE Read_Options(FileList, UseCol, InfoFileName, DistMatrixName, &
   NumArg, c, simplify, UseCut, k_value, AnalyticGrad, Plot, Saveplotdata, Step, &
   ConvCriterion, MaxIterations, MaxNumConf,RdBigMat,WrBigMat, BigMatFile)
!*****************************************************************************************
!This subroutine reads all the options that can be given to the
!program through the command line
!-----------------------------------------------------------------------------------------
   use parameters
   implicit none
   !-VARIABLE DEFINITIONS--------------------------------------------------------------
   integer MaxNumConf                                      !Max number of conf. states by torsion allowed
   logical GivenInfoFileName                               !The info filename is given? (True/False)
   logical GivenDistMatrixName                             !The distance matrix filename is given? (T/F)
   logical UseCut                                          !Use cutoff criterion? (True/False)
   logical Plot                                            !Plot data (True/False)
   logical Saveplotdata                                    !Save Plot data (True/False)
   character*60 DistMatrixName                             !Distance matrix filename
   character*60 InfoFileName                               !Info-file filename
   integer NumArg                                          !Number of arguments in the command line
   character*60 Arguments(NumArg)                          !Arguments in the command line
   character*60 FileList(NumArg)                           !List of data files (*.dat)
   character*3 simplify                                    !Eliminate rigid torsions? (True/False)
   real(DP) k_value                                        !k_value (see default values in the main program)
   character*3 answer                                      !Use analytic gradient? (yes/no)
   logical  AnalyticGrad                                   !The logical analogous to "answer" (True/False)
   real(DP) Step                                           !Step size in the optimization
   integer  MaxIterations                                  !Max number of iterations in the optimization
   real(DP) ConvCriterion                                  !Convergence criterion for the optimization
   integer UseCol                                          !Column to be used for the discretization in the
   !data files
   logical RdBigMat,WrBigMat
   character*120  BigMatFile

   real(DP)  RealVar                                       !
   integer i, c                                             !
   integer ios !                                            >  Auxiliary or Dummy variables
   logical EO                                              !
   character*60 arg                                        !
   !-----------------------------------------------------------------------------------


!-----INITIAL VALUES
   GivenInfoFileName=.false.
   GivenDistMatrixName=.false.
   EO=.false.                                              !End Options is False

!--READING OPTIONS
   do i = 1, NumArg
      call getarg(i, arg)
      Arguments(i)=arg
   enddo

   i = 1
   c = 0
   do while (i .le. NumArg)
      if((Arguments(i).eq.'-u').or.(Arguments(i).eq.'-usecol').and.(.not.EO)) then
         read(Arguments(i+1), *,iostat = ios) UseCol
         read(Arguments(i+1), *,iostat = ios) RealVar
         if((ios .ne. 0).or.(UseCol .lt. 1).or.(dfloat(UseCol) .ne. RealVar).and.(.not.EO)) &
            stop "ERROR: Check the-u/-usecol option. Use-help option for quick help"
         i = i+1
      elseif((Arguments(i).eq.'-dist').and.(.not.EO)) then
         read(Arguments(i+1), '(A60)',iostat = ios) DistMatrixName
         if(ios .ne. 0) &
            stop "ERROR: Check the file name of distance matrix. Use-help option for quick help"
         i = i+1
         GivenDistMatrixName=.true.
      elseif((Arguments(i).eq.'-i').or.(Arguments(i).eq.'-info').and.(.not.EO)) then
         read(Arguments(i+1), '(A60)',iostat = ios) InfoFileName
         if(ios .ne. 0) stop "ERROR: Check the info file name. Use-help option for quick help"
         i = i+1
         GivenInfoFileName=.true.
      elseif((Arguments(i).eq.'-nocut').and.(.not.EO)) then
         UseCut=.false.
      elseif((Arguments(i).eq.'-s').or.(Arguments(i).eq.'-simplify').and.(.not.EO)) then
         read(Arguments(i+1), *,iostat = ios) simplify
         if((ios .ne. 0).or.(.not.((simplify.eq."yes").or.(simplify.eq."no")))) &
            stop "ERROR: Check the option simplyfy. Use-help option for quick help"
         i = i+1
         GivenInfoFileName=.true.
      elseif((Arguments(i).eq.'-k').and.(.not.EO)) then
         call getarg(i+1, arg)
         read(arg, *,iostat = ios) k_value
         if(ios > 0) stop "ERROR: Check k_value. Use-help option for quick help"
         if(k_value .le. 0) stop "ERROR: k_value most be a positive real"
         i = i+1
      elseif((Arguments(i).eq.'-wrbigmat').and.(.not.EO)) then
         read(Arguments(i+1), *,iostat = ios) BigMatFile
         if(ios .ne. 0) then
            print*,"ERROR: Check the option -wrbigmat "
            print*,"Use-help option for quick help"
            stop
         endif
         WrBigMat=.true.
         i = i+1
      elseif((Arguments(i).eq.'-rdbigmat').and.(.not.EO)) then
         read(Arguments(i+1), *,iostat = ios) BigMatFile
         if(ios .ne. 0) then
            print*,"ERROR: Check the option -rdbigmat "
            print*,"Use-help option for quick help"
            stop
         endif
         RdBigMat=.true.
         i = i+1
      elseif((Arguments(i).eq.'-ag').and.(.not.EO)) then
         read(Arguments(i+1), *,iostat = ios) answer
         if((ios .ne. 0).or.(.not.((answer.eq."yes").or.(answer.eq."no")))) then
            print*,"ERROR: Check the option-ag. "
            print*,"Use-help option for quick help"
            stop
         endif
         if(answer.eq."yes") AnalyticGrad=.true.
         i = i+1
      elseif((Arguments(i).eq.'-saveplotdata').and.(.not.EO)) then
         read(Arguments(i+1), *,iostat = ios) answer
         if((ios .ne. 0).or.(.not.((answer.eq."yes").or.(answer.eq."no")))) then
            print*,"ERROR: Check the option-saveplotdata "
            print*,"Use-help option for quick help"
            stop
         endif
         if(answer.eq."yes") Saveplotdata=.true.
         i = i+1
      elseif((Arguments(i).eq.'-plot').and.(.not.EO)) then
         read(Arguments(i+1), *,iostat = ios) answer
         if((ios .ne. 0).or.(.not.((answer.eq."yes").or.(answer.eq."no")))) then
            print*,"ERROR: Check the option-plot "
            print*,"Use-help option for quick help"
            stop
         endif
         if(answer.eq."yes") Plot=.true.
         i = i+1
      elseif((Arguments(i).eq.'-step').and.(.not.EO)) then
         call getarg(i+1, arg)
         read(arg, *,iostat = ios) Step
         if(ios > 0) stop "ERROR: Check Step value. Use-help option for quick help"
         if(Step .le. 0.0d0) stop "ERROR: Step most be a positive real"
         i = i+1
      elseif((Arguments(i).eq.'-crit').and.(.not.EO)) then
         call getarg(i+1, arg)
         read(arg, *,iostat = ios) ConvCriterion
         if(ios > 0) stop "ERROR: Check Convergence Criterion value. Use-help option for quick help"
         if(ConvCriterion .le. 0) stop "ERROR: Convergence Criterion most be a positive real"
         i = i+1
      elseif((Arguments(i).eq.'-maxitr').and.(.not.EO)) then
         read(Arguments(i+1), *,iostat = ios) MaxIterations
         read(Arguments(i+1), *,iostat = ios) RealVar
         if((ios .ne. 0).or.(MaxIterations .lt. 1).or.(MaxIterations .ne. RealVar)) &
            stop "ERROR: Check the-maxitr option. Use-help option for quick help"
         i = i+1
      elseif((Arguments(i).eq.'-maxconf').and.(.not.EO)) then
         read(Arguments(i+1), *,iostat = ios) MaxNumConf
         read(Arguments(i+1), *,iostat = ios) RealVar
         if((ios .ne. 0).or.(MaxNumConf .lt. 2).or.&
            (MaxNumConf .ne. RealVar).or.(MaxNumConf .gt. 9)) then
            print*,"ERROR: Check the-maxconf option. "
            print*,"The maximum number of conformers to select from"
            print*,"each dihedral, most be an integer number between 2 and 9."
            print*,"Use-help option for quick help"
            stop
         endif
         i = i+1
      elseif((Arguments(i).eq.'-help').or.(Arguments(i).eq.'-h').and.(.not.EO)) then
         print*,"QUICK HELP:"
         print*,""
         print*,"USAGE: centro_prep [OPTIONS] file1 file2 ..."
         print*," "
         print*,"The files file1, file2 etc, contain the time series of the dihedral angles."
         print*,"It is also valid the usage of regular expressions, for example:"
         print*,"             "
         print*," centro_prep [OPTIONS] file?"
         print*," centro_prep [OPTIONS] fi*"
         print*,"             "
         print*,"OPTIONS:"
         print*,"--------"
         print*," -u/-usecol NUMCOL                             Default:2"
         print*,"            The number NUMCOL specifies which column of file1, file2 etc."
         print*,"            contains the dihedral angle variable. This value is normally 2"
         print*,"            since the first column is often the time or the snapshot number"
         print*,"      "
         print*," -dist      DISTANCE_MATRIX_FILE_NAME          Default:distance_matrix.dat"
         print*,"            Specifies the distance matrix file name of the full system in"
         print*,"            vacuo."
         print*,""
         print*," -i/-info   DIH_INFO_FILE_NAME                 Default:atoms_in_tor.info"
         print*,"            Specifies the file name that contains which atoms are involved"
         print*,"            in each dihedral (the two central ones). For example, if the "
         print*,"            first row is 3 4, means that the first dihedral (the one whose"
         print*,"            time series is in file1) is defined over the bond 3-4."
         print*,""
         print*," -nocut                                        Default:Use cut-off"
         print*,"            Using this option no cut-off will be applied and the options"
         print*,"            (-dist/-info) and files DISTANCE_MATRIX_FILE_NAME and"
         print*,"            DIH_INFO_FILE_NAME  are not needed."
         print*,""
         print*," -s/-simplify  yes/no                          Default:yes"
         print*,"            Removes the quasifrozen discrete dihedrals with hardly any "
         print*,"            conformational change."
         print*,""
         print*," -k         K_VALUE                            Default: 0.5"
         print*,"            The k_value sets the smoothing parameter v "
         print*,"            in the von-Mises kernel density estimation as"
         print*,"            proposed in eq.(7) of ref: "
         print*,"            Computational Statistics & Data Analysis"
         print*,"            Volume 52, Issue 7, 15 March 2008, Pages 3493-3500."
         print*,"            We set by default k_value = 0.5 because it slightly"
         print*,"            oversmooths the PDFs, which is convenient for minimizations."
         print*,""
         print*," -ag        yes/no                             Default: no"
         print*,"            Use analytic gradient on each step of the minimizations"
         print*,"            that search minima in the Probability Density"
         print*,"            Functions (PDFs). The default option (no) uses a quite"
         print*,"            accurate and fast linear interpolation of the gradient."
         print*,""
         print*," -step      STEP_SIZE                          Default: 5 (degrees)"
         print*,"            Step size a minimization of the form:"
         print*,"            X(n+1)=X(n)-STEP_SIZE*GRADIENT"
         print*,""
         print*," -crit      CONVERGENCE_CRITERION              Default: 1.0E-4"
         print*,"            Convergence criterion for the minimization"
         print*,""
         print*," -maxitr    MAX_NUMBER_OF_ITERATIONS           Default: 1000"
         print*,"            Maximum number of iteration in each minimization"
         print*,""
         print*," -maxconf   MAX_NUMBER_OF_CONFORMERS           Default: 3"
         print*,"            Maximum number of conformers to be searched for each"
         print*,"            dihedral angle"
         print*,""
         print*," -plot      yes/no                             Default: no"
         print*,"            Plots PDF and time evolution of Dihedral angles "
         print*,"            The DISLIN package is used to create PNG files"
         print*,""
         print*," -saveplotdata   yes/no                       Default: no"
         print*,"            Plot data are saved in csv files."
         print*,"            Useful to recreate graphics using other software."
         print*,""
         print*," --wrbigmat BIGMATRIX_FILE_NAME                Default: no "
         print*,"            The BIGMATRIX is printed out and the program stops."
         print*,"            This option can be useful to distribute the processing  "
         print*,"            of many (and large) dihedral files across several nodes/proc."
         print*,""
         print*," --rdbigmat BIGMATRIX_FILE_NAME                Default: no "
         print*,"            The BIGMATRIX is read (not built) and then the program continues"
         print*,"            doing satistical analysis and constructing reduced_dist_matrix.dat "
         print*,"            This option can be useful to distribute the processing  "
         print*,"            of many (and large) dihedral files across several nodes/proc."
         print*,""
         print*," -help"
         print*,"            Prints this quick help"
         print*,""
         print*,"EXAMPLE:"
         print*,"---------"
         print*,"cencalc_prep  -ag yes -k 0.5 -plot.yes d????.dat"
         stop
      else
         EO=.true.
         arg = Arguments(i)
         if(arg(1:1).eq."-") &
            stop "ERROR: Check the options, use-help option for quick help"
         c = c+1
         FileList(c)=arg
      endif
      i = i+1
   enddo

   return
END SUBROUTINE Read_Options

!*****************************************************************************************
SUBROUTINE HistPDF_Time_Plot(Saveplotdata, NumSnap, dihedFILE, dihedID, FNumMin, minimum, &
   Pop, DiscreteAng, density, histogram, plotsize)
!*****************************************************************************************
! Get the plots of the Von Mises probability density function
! The DISLIN package is needed
!
!  Things yet to do: Implement control of quality PNG
!----------------------------------------------------------------------------------------
   use parameters
   use dislin
   implicit none
   !-VARIABLE DEFINITIONS--------------------------------------------------------------
   logical  Saveplotdata          ! Save plot data in a csv file
   integer, intent(in):: NumSnap                          !Number of snapshots
   integer, intent(in):: DiscreteAng(NumSnap)      !Discretized Dihedral angle
   integer plotsize               ! Scaling factor for plot size (1, 2, 3, ....
   integer  FnumMin               ! Number of located minima
   real(DP) minimum(9)            ! Positining of minima (in degrees)
   real(DP) Pop(9)                ! % population
   real(DP) density(0:360)        ! von Mises Probability Density Function
   real(DP) histogram(0:72)       ! Histogram Probability Density Function 5 degree binsize
   character*60  dihedFILE        ! dihedral angle Filename
   character*60  dihedID          ! ID of dihedral angle
!
   real(DP)  xsnap(0:1980)  ! For plotting time evolution
   real(DP)  ysnap(0:1980)  ! For plotting time evolution
   real(DP) Ang(0:360), AngBin(0:72), Y1Bin(0:72)
   character legtext*60, legmin*60
   integer  i, ic, lID, lFILE, ix, istep

!     Pruning the DiscreteAng array for plotting
   do i = 1, 1980
      xsnap(i)=0.0d0
      ysnap(i)=0.0d0
   enddo
   if ( NumSnap .le. 1980) then
      do i = 1, NumSnap
         xsnap(i) = dfloat(i)
         ysnap(i) = dfloat(DiscreteAng(i))
      enddo
      ix = NumSnap
   else
      istep = INT(NumSnap/1980)+1
      ix = 0
      do  i = 1, NumSnap, istep
         ix = ix+1
         xsnap(ix)=dfloat(i)
         ysnap(ix)=dfloat(DiscreteAng(i))
      enddo
   endif
!     Setting auxiliary arrays for PDF plots
   do i = 0, 360
      Ang(i)=dfloat(i)
   enddo
   do i = 0, 72
      AngBin(i)=dfloat(i)*5.0d0
      Y1Bin(i)=0.0d0
   enddo
!     Finding the actual length of dihedFILE
   i = 1
   do while ((dihedFILE(i:i) .ne. '') .and. (dihedFILE(i:i) .ne. ' ') .and. (i .le. 60))
      i = i+1
   enddo
   lFILE = i-1

!
! General settings of the DISLIN plot
!
   CALL METAFL('PNG')      !   PNG file type
   CALL WINSIZ(853*plotsize, 603*plotsize)    !  resolution of image file
!     CALL PNGMOD('ON','TRANSPARENCY')
   CALL FILMOD('DELETE')
   CALL SETFIL(dihedFILE(1:lFILE)//'.png')   ! PNG file name

!     CALL METAFL('TIFF')      !   file type
!     CALL WINSIZ(853*plotsize, 603*plotsize)    !  resolution of image file
!     CALL FILMOD('DELETE')
!     CALL SETFIL(dihedFILE(1:lFILE)//'.tif')

   CALL SCRMOD('REVERS')
   CALL SETPAG('DA4L')    ! This is default A4  2970 x 2100

   CALL DISINI
   CALL PAGERA            !  Prints the page frame
   CALL COMPLX            !  Sets a complex (nicer) font
   CALL BMPFNT('COMPLEX')
   CALL SOLID             ! solid lines
   CALL LINWID(plotsize*1)  ! line width
   CALL SCLFAC(dfloat(plotsize)*1.0d0)
   CALL SCLMOD('FULL')

   CALL TITLIN  &
      ('PDF and Time Plots for '//dihedFILE(1:lFILE)//'='//dihedID, 1)
   CALL HEIGHT(35)           ! Height of characters used for Title and labels
   CALL LINESP(0.10d0)       ! 0.10*HEIGHT defines the spacing between Title and Axis
!
!    Setting the PDF plot
!
   CALL HNAME(36)           ! Height of characters used for Axis names
   CALL NAME('Dihed Angle (degree)','X')
   CALL NAME('Unnormalized PDF','Y')
   CALL LABDIG(-1, 'X')    !  None decimal places for X labels
   CALL LABDIG( 1, 'Y')    !  1 decimal places for Y labels
   CALL AXSPOS(450, 1000)  !  Positioning of the left down corner
   CALL AXSLEN(1980, 800)  !  Lenght of Graf Axes
   CALL AXSBGD(-1)        !  Background color of Graf Area NONE
!
!     Printing the label for the legend box (Minima positioning)
   write(legmin, '(6(1X, F5.1))') (minimum(i), i = 1, FNumMin)
   CALL LEGINI(legtext, 1, 30)
   CALL LEGTIT('')
   CALL LEGLIN(legtext, 'Minima: '//legmin, 1)
!
!     GRAF (X0, X1, X of First label, Label space, Idem for Y axis)
!
   CALL GRAF(0.D0, 360.D0, 0.D0, 30.D0, 0.D0, 1.5D0, 0.D0, 0.25D0)
   CALL SETRGB(0.7D0, 0.7D0, 0.7D0)
   CALL GRID(1, 1)           !  Number of grid lines between labels

   CALL COLOR('FORE')  ! Resets color to default value
   CALL TITLE         ! Print Title  OVER the Axis System

   CALL NOCHEK

   CALL COLOR('RED')
   CALL THKCRV(2)
   CALL CURVE(Ang, density, 361)

   CALL LINESP(0.50d0)
   CALL COLOR('FORE')
   CALL LEGEND(legtext, 8)

   CALL COLOR('BLUE')
   CALL BARS(AngBin, Y1Bin, histogram, 72)

   CALL ENDGRF
!
!     Setting the time evolution plot
!
   CALL COLOR('FORE')  ! Resets color to default value
   CALL HNAME(36)           ! Height of characters used for Axis names
   CALL NAME('NumSnap','X')
   CALL NAME('Conf. State.','Y')
   CALL LABDIG(-1, 'X')    !  None decimal places for X labels
   CALL LABDIG(-1, 'Y')    !  1 decimal places for Y labels
   CALL AXSPOS(450, 1900)  !  Positioning of the left down corner
   CALL AXSLEN(1980, 700)  !  Lenght of Graf Axes
   CALL AXSBGD(-1)        !  Background color of Graf Area NONE

!     Printing the label for the legend box (Minima positioning)
   write(legmin, '(6(1X, F5.1))') (Pop(i), i = 1, FNumMin)
   CALL LEGINI(legtext, 1, 30)
   CALL LEGTIT('')
   CALL LEGLIN(legtext, '% Pop.: '//legmin, 1)
!
   CALL GRAF(0.D0,dfloat(NumSnap), 0.D0, &
      dfloat(NumSnap)/5.0d0, 0.D0, 4.D0, 0.D0, 1.D0)
   CALL SETRGB(0.7D0, 0.7D0, 0.7D0)
   CALL GRID(1, 0)           !  Number of grid lines between labels

   CALL NOCHEK

   CALL COLOR('GRAY')
   CALL THKCRV(1)
   CALL CURVE(xsnap, ysnap, ix)
!
   CALL LINESP(0.50d0)
   CALL COLOR('FORE')
   CALL LEGEND(legtext, 8)

   CALL ENDGRF

   CALL ERRMOD('protocol','off')

   CALL DISFIN
   if ( Saveplotdata ) then
      open(33, file = dihedFILE(1:lFILE)//'.csv',status='unknown')
      write(33, '(a)') '# Dihedral= '//dihedID
      write(33, '(''# Density plot (0:360) '')')
      write(33, '(''# Ang, Density '')')
      do i = 0, 360
         write(33, '(F10.5, '', '',F10.5)') Ang(i), density(i)
      enddo
      write(33, '(''# Hist plot (0:72) '')')
      write(33, '(''# Angbin, Histogram'')')
      do i = 0, 72
         write(33, '(F10.5, '', '',F10.5)') Angbin(i), histogram(i)
      enddo
      write(33, '(''# Time evolution plot '')')
      write(33, '(''# X, Y '')')
      do i = 1, ix
         write(33, '(F20.1, '', '',F5.1)') xsnap(i), ysnap(i)
      enddo
      close(33)
   endif
!
   return
END SUBROUTINE HistPDF_Time_Plot


!*****************************************************************************************
SUBROUTINE Get_PDF_Min(ang, NumSnap, k_value, &
   AnalyticGrad, Step, ConvCriterion, MaxIterations, &
   MaxNumConf, minimum, FNumMins, density, histogram)
!*****************************************************************************************
!  Builds Histogram and von-Mises PDF and determines the positioning of the
!  PDF minima
!----------------------------------------------------------------------------------------
   use parameters
   implicit none
   !-VARIABLE DEFINITIONS--------------------------------------------------------------
   integer, intent(in):: NumSnap                          !Number of snapshots
   real(DP), intent(in):: ang(NumSnap)                     !Dihedral angle
   real(DP), intent(in):: k_value                          !k_value (see default values in the main program)
   logical, intent(in):: AnalyticGrad                     !Use analytic gradient? (true/false)
   real(DP), intent(in):: Step                             !Step for the steepest-descendent optimization
   real(DP), intent(in):: ConvCriterion                    !Convergence criterion for the steepest-descendent opt.
   integer, intent(in):: MaxIterations                    !Max number of iterations in the optimization
   integer, intent(in):: MaxNumConf                       !Max number of conformers allowed by torsion
   real(DP), intent(out):: minimum(9)                      !Position of the minimums in the PDF
   integer, intent(out):: FNumMins                        !Final number of minimums founded
   real(DP), intent(out):: density(0:360)                  ! Von Mises Probability Density Function (PDF)
   real(DP), intent(out):: histogram(0:72)                 ! Histogram Probability Density Function (PDF) 5 degree binsize
   real(DP) CoordMaxs(9)                                   !Positions of the maximums in the PDF
   real(DP) SecDerivMaxs(9)                                !Positions of the maximums of the 2nd derivate of the (PDF)
   real(DP) BessI                                          !Modified Bessel function
   real(DP) gradient(0:360)                                !Gradient of the PDF
   real(DP) secderiv(0:360)                                !Second derivate of the PDF
   real(DP) SmoothingParam                                 !Smoothing parameter (see ref. 2)
   real(DP) vMgradient                                     !Gradient in a general position (not only in {0, ...359})
   real(DP) vMdensity                                      !Density in a general position
   real(DP) vMsec_deriv                                    !Second derivate in a general position
   real(DP) Func_vMgradient                                !Function
   real(DP) Step_times_grad                                !Step x gradient
   integer ObsNmax                                         !Observed Number of Maximums

   integer contador                                        !
   integer i, j, Itrs, m                                      !
   real(DP) temp                                           !> Auxiliary or dummy variables
   logical Converged                                       !
   real(DP) grad                                           !
!-----------------------------------------------------------------------------------
!-----Computing Histogram
   do j = 0, 72
      histogram(j)=0.0d0
   enddo
   do j = 0, 360
      density(j)=0.0d0
   enddo
   do i = 1, NumSnap
      j = INT(ang(i)/5.0d0)
      histogram(j)=histogram(j)+1.0d0
   enddo
   do j = 0, 72
      histogram(j)= (5.0d0)*histogram(j)/dfloat(NumSnap)
   enddo
!
   secDerivMaxs = 0.d0
!----Computing the Bandwith
   SmoothingParam=(3.d0*dfloat(NumSnap)*(k_value**2)*BessI(2, 2*k_value)/&
      (4.0d0*(Pi**0.5d0)*(BessI(0, k_value))**2.d0))**(2.d0/5.d0)


!----Calculating densities, gradient, second derivate
!----and looking for local maximums
   ObsNmax = 0    ! Observed Number of Maximums
   do i = 0, 360
      call vMises(NumSnap, SmoothingParam, (dfloat(i)*Pi/180.d0), ang, vMdensity, vMgradient, vMsec_deriv)
      gradient(i)=vMgradient
      density(i)=vMdensity
      secderiv(i)=vMsec_deriv
      if(i .gt. 0) then
         if(((gradient(i-1)*gradient(i)).le.0.d0).and.(vMsec_deriv.lt.(-1.0D-4)).and.(ObsNmax .lt. 9)) then
            ObsNmax = ObsNmax+1
            CoordMaxs(ObsNmax)=dfloat(i)/2.d0+dfloat(i-1)/2.d0           !Approximate coord of the maximum
            SecDerivMaxs(ObsNmax)=&
               secderiv(i)/2.d0 + secderiv(i-1)/2.d0              !Approximate second derivate
         endif
      endif
   enddo

!----Sorting to give preference to those maxima with
!----lower(more negative) second derivates
   do i = 1, ObsNmax-1
      do j = i+1, ObsNmax
         if(SecDerivMaxs(i).gt.SecDerivMaxs(j)) then
            temp = SecDerivMaxs(i)
            SecDerivMaxs(i)=SecDerivMaxs(j)
            SecDerivMaxs(j)=temp
            temp = CoordMaxs(i)
            CoordMaxs(i)=CoordMaxs(j)
            CoordMaxs(j)=temp
         endif
      enddo
   enddo

!----In a circular variable, the Max. number of minima
!----most be equal to the Max. number of maxima, we get
!----the lowest value between the predefined MaxNumConf
!----and the observed number of maxima/minima (ObsNmax)
   if(MaxNumConf .lt. ObsNmax)  FNumMins = MaxNumConf
   if(MaxNumConf .ge. ObsNmax)  FNumMins = ObsNmax

!----Sorting only the first FNumMins using the coordinates
   do i = 1, FNumMins-1
      do j = i+1, FNumMins
         if(CoordMaxs(i).gt.CoordMaxs(j)) then
            temp = CoordMaxs(i)
            CoordMaxs(i)=CoordMaxs(j)
            CoordMaxs(j)=temp
         endif
      enddo
   enddo

!----Guessing the minima locations
   if (FNumMins .eq. 1) then
      minimum(1)=(CoordMaxs(1)-180.0d0)
      if (minimum(1).lt.0.0d0) minimum(1)=minimum(1)+360.d0
   else
      do i = 1, FNumMins
         if (i .eq. 1) then
            minimum(i)=(CoordMaxs(i)+CoordMaxs(FNumMins)-360.d0)/2.d0
            if (minimum(i).lt.0.0d0) minimum(i)=minimum(i)+360.d0
         else
            minimum(i)=(CoordMaxs(i)+CoordMaxs(i-1))/2.d0
         endif
      enddo
   endif

!---Optimizing Minima
   if(FNumMins .gt. 1) then
      do m = 1, FNumMins
         Itrs = 0
         Converged=.false.
         do while((.not.Converged).and.(Itrs < MaxIterations))
            Itrs = Itrs+1
            if(AnalyticGrad) then
               grad = Func_vMgradient(NumSnap, SmoothingParam, (Pi/180.d0)*minimum(m), ang)
            else
               grad=(abs(minimum(m)-1.d0-dfloat( int(minimum(m))) )   )&        !Gradient Linear interpolation
                  *gradient(int(minimum(m)))+&
                  (abs(minimum(m)-  dfloat(int(minimum(m)))   )&
                  *gradient(1+int(minimum(m))))
            endif
            if ((abs(grad).lt.ConvCriterion)) then
               Converged=.true.
            else
               if(abs(Step*grad).gt.5.0d0) then
                  Step_times_grad=((Step*grad)/abs(Step*grad))*5.0d0
               else
                  Step_times_grad = Step*grad
               endif
               minimum(m)=minimum(m)-Step_times_grad
               if(minimum(m).lt.0.d0) minimum(m)=minimum(m)+360.d0
               if(minimum(m).ge.360.d0) minimum(m)=minimum(m)-360.d0
            endif
         enddo
      enddo
   endif

   do i = 1, FNumMins-1
      do j = i+1, FNumMins
         if(minimum(i).gt.minimum(j)) then
            temp = minimum(i)
            minimum(i)=minimum(j)
            minimum(j)=temp
         endif
      enddo
   enddo


   return


END SUBROUTINE Get_PDF_Min
!----------------------------------------------------------------------------------------



!*****************************************************************************************
SUBROUTINE vMises(NumSnap, SmoothingParam, X, ang, vMdensity, vMgradient, vMsec_deriv)
!*****************************************************************************************
!Computes the density, gradient and second derivate in a position X
!-----------------------------------------------------------------------------------------
   use parameters
   implicit none
   !-VARIABLE DEFINITIONS--------------------------------------------------------------
   integer NumSnap                                         !Number of snapshots
   real(DP) vMdensity                                      !von Mises density
   real(DP) vMgradient                                     !von Mises gradient
   real(DP) vMsec_deriv                                    !von Mises second derivate
   real(DP) BessI                                          !Modified Bessel function
   real(DP) SmoothingParam                                 !Smoothing parameter
   real(DP) ang(*)                                         !Data Angles in degrees
   real(DP) aj                                             !Data Angles in radians
   real(DP) expo                                           !Auxiliary real variable
   integer j                                               !Dummy integer variable
   real(DP), intent(in):: X                               !The position where we want to
   !evaluate vMdensity, vMgradient, ...
   !-----------------------------------------------------------------------------------

   vMgradient = 0.0d0
   vMgradient = 0.0d0
   vMsec_deriv = 0.0d0
   do j = 1, NumSnap
      aj = ang(j)*(Pi/180.d0)                                    !transforming to radians
      expo = exp( SmoothingParam*cos(X-aj))
      vMdensity = vMdensity+expo
      vMgradient = vMgradient-sin(X-aj)*expo
      vMsec_deriv = vMsec_deriv+&
         (SmoothingParam*(sin(X-aj))**2.d0-cos(X-aj))*expo
   enddo
   vMdensity = vMdensity/( 2.0d0*Pi*dfloat(NumSnap)*BessI(0, SmoothingParam))
   vMgradient = vMgradient*SmoothingParam/( 2.0d0*Pi*dfloat(NumSnap)*BessI(0, SmoothingParam))
   vMsec_deriv = vMsec_deriv*SmoothingParam/( 2.0d0*Pi*dfloat(NumSnap)*BessI(0, SmoothingParam))
   return
END SUBROUTINE vMises
!----------------------------------------------------------------------------------------



!*****************************************************************************************
FUNCTION Func_vMgradient(NumSnap, SmoothingParam, X, ang)
!*****************************************************************************************
!Function specifically designed to compute the von Misses Gradient
!----------------------------------------------------------------------------------------
   use parameters
   implicit none
   !-VARIABLE DEFINITIONS--------------------------------------------------------------
   real(DP) Func_vMgradient                                !
   real(DP) BessI                                          !Modified Bessel function
   real(DP) SmoothingParam                                 !Smoothing parameter
   real(DP) ang(*)                                         !Data Angle in degrees
   real(DP) aj                                             !Data Angle in radians
   integer NumSnap                                         !Number of snapshots
   integer j                                               !
   real(DP), intent(in):: X                               !The position where we want to
   !evaluate vMdensity, vMgradient, ...
   !-----------------------------------------------------------------------------------
   Func_vMgradient = 0.0d0
   do j = 1, NumSnap
      aj = ang(j)*(Pi/180.d0)                                    !transforming to radians
      Func_vMgradient = Func_vMgradient-SmoothingParam*sin(X-aj)*&
         exp( SmoothingParam*cos(X-aj))/&
         ( 2.0d0*Pi*dfloat(NumSnap)*BessI(0, SmoothingParam))
   enddo
   return
END FUNCTION
!----------------------------------------------------------------------------------------



!*****************************************************************************************
FUNCTION BessI(N, X)
!*****************************************************************************************
!     Calc. the first kind modified Bessel function
!
!     REFERENCE:
!     C.W.CLENSHAW, CHEBYSHEV SERIES FOR MATHEMATICAL FUNCTIONS,
!     MATHEMATICAL TABLES, VOL.5, 1962.
!-----------------------------------------------------------------------------------------
   use parameters
   implicit none
   integer IACC, N, I, J , M
   real(DP) X, BessI, BIGNO, BIGNI, BessI0, BessI1, TOX, BIM, BI, BIP
   parameter (IACC = 40, BIGNO = 1.D10, BIGNI = 1.D-10)

   if (N .EQ. 0) then
      BessI = BessI0(X)
      return
   endif

   if (N .EQ. 1) then
      BessI = BessI1(X)
      return
   endif

   if(X .EQ. 0.D0) then
      BessI = 0.D0
      return
   endif

   TOX = 2.D0/X
   BIP = 0.D0
   BI  = 1.D0
   BessI = 0.D0
   M = 2*((N+INT(SQRT(DFLOAT(IACC*N)))))

   do J = M, 1, -1
      BIM = BIP+DFLOAT(J)*TOX*BI
      BIP = BI
      BI  = BIM
      if (ABS(BI).GT.BIGNO) then
         BI  = BI*BIGNI
         BIP = BIP*BIGNI
         BessI = BessI*BIGNI
      endif
      if (J .eq. N) BessI = BIP
   enddo

   BessI = BessI*BessI0(X)/BI

   return

END FUNCTION BessI
!----------------------------------------------------------------------------------------



!*****************************************************************************************
FUNCTION BessI0(X)
!*****************************************************************************************
! Auxiliary Bessel functions for N = 0
!-----------------------------------------------------------------------------------------
   use parameters
   real(DP) X, BessI0, Y, P1, P2, P3, P4, P5, P6, P7,  &
      Q1, Q2, Q3, Q4, Q5, Q6, Q7, Q8, Q9, AX, BX
   data P1, P2, P3, P4, P5, P6, P7/1.D0, 3.5156229D0, 3.0899424D0, 1.2067429D0,  &
      0.2659732D0, 0.360768D-1, 0.45813D-2/
   data Q1, Q2, Q3, Q4, Q5, Q6, Q7, Q8, Q9/0.39894228D0, 0.1328592D-1, &
      0.225319D-2, -0.157565D-2, 0.916281D-2, -0.2057706D-1,  &
      0.2635537D-1, -0.1647633D-1, 0.392377D-2/

   if(ABS(X).LT.3.75D0) then
      Y=(X/3.75D0)**2
      BessI0 = P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7)))))
   else
      AX = ABS(X)
      Y = 3.75D0/AX
      BX = EXP(AX)/SQRT(AX)
      AX = Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))))
      BessI0 = AX*BX
   endif

   return
END FUNCTION BessI0
!-----------------------------------------------------------------------------------------



!*****************************************************************************************
FUNCTION BessI1(X)
!*****************************************************************************************
! Auxiliary Bessel functions for N = 1
!-----------------------------------------------------------------------------------------
   use parameters
   real(DP) X, BessI1, Y, P1, P2, P3, P4, P5, P6, P7,  &
      Q1, Q2, Q3, Q4, Q5, Q6, Q7, Q8, Q9, AX, BX
   data P1, P2, P3, P4, P5, P6, P7/0.5D0, 0.87890594D0, 0.51498869D0,  &
      0.15084934D0, 0.2658733D-1, 0.301532D-2, 0.32411D-3/
   data Q1, Q2, Q3, Q4, Q5, Q6, Q7, Q8, Q9/0.39894228D0, -0.3988024D-1, &
      -0.362018D-2, 0.163801D-2, -0.1031555D-1, 0.2282967D-1, &
      -0.2895312D-1, 0.1787654D-1, -0.420059D-2/

   if(ABS(X).LT.3.75D0) then
      Y=(X/3.75D0)**2
      BessI1 = X*(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
   else
      AX = ABS(X)
      Y = 3.75D0/AX
      BX = EXP(AX)/SQRT(AX)
      AX = Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))))
      BessI1 = AX*BX
   endif

   return
END FUNCTION BessI1
!-----------------------------------------------------------------------------------------


!*****************************************************************************************
SUBROUTINE SORT2(N, RA, IB)
!*****************************************************************************************
!     SORTS AN ARRAY RA OF LENGTH N INTO ASCENDING NUMERICAL ORDER
!     USING THE HEAPSORT ALGORITHM, WHILE MAKING THE CORRESPONDING
!     REARRANGEMENT OF THE ARRAY IB.
!
   use parameters
   implicit none
   !-VARIABLE DEFINITIONS--------------------------------------------------------------
   real(DP) RA(*)
   integer  N, IB(*)
   integer I, J, L, IR, IRB
   real(DP) RRA, RRB
   IF ( ( N .EQ. 0 ) .OR. ( N .EQ. 1 ) )  RETURN
   L = N/2+1
   IR = N
10 CONTINUE
   IF(L .GT. 1)THEN
      L = L-1
      RRA = RA(L)
      IRB = IB(L)
   ELSE
      RRA = RA(IR)
      IRB = IB(IR)
      RA(IR)=RA(1)
      IB(IR)=IB(1)
      IR = IR-1
      IF(IR .EQ. 1)THEN
         RA(1)=RRA
         IB(1)=IRB
         RETURN
      ENDIF
   ENDIF
   I = L
   J = L+L
20 IF(J .LE. IR)THEN
      IF(J .LT. IR)THEN
         IF(RA(J).LT.RA(J+1))J = J+1
      ENDIF
      IF(RRA .LT. RA(J))THEN
         RA(I)=RA(J)
         IB(I)=IB(J)
         I = J
         J = J+J
      ELSE
         J = IR+1
      ENDIF
      GO TO 20
   ENDIF
   RA(I)=RRA
   IB(I)=IRB
   GO TO 10
   return
END SUBROUTINE SORT2
!
SUBROUTINE COLSPEC(N,INDX,L,KEY)
!*****************************************************************************************
!
   use parameters
   implicit none
   !-VARIABLE DEFINITIONS--------------------------------------------------------------
   integer N
   integer INDX(N)
   character KEY*512,TMPKEY*512
   character TERM*5
   integer I,J,L
   logical consecutive

   WRITE(KEY,'(I4)') INDX(1)
   L=4
   DO I=2,N
      IF (INDX(I)-INDX(I-1) .gt. 1) THEN
         IF ( I .eq. 2 ) THEN
            write(TERM,'('','',I4)') INDX(I)
            KEY=KEY(1:L)//TERM
            L=L+5
         ELSE
            IF ( consecutive ) THEN
               write(TERM,'(''-'',I4)') INDX(I-1)
               KEY=KEY(1:L)//TERM
               L=L+5
               write(TERM,'('','',I4)') INDX(I)
               KEY=KEY(1:L)//TERM
               L=L+5
            ELSE
               write(TERM,'('','',I4)') INDX(I)
               KEY=KEY(1:L)//TERM
               L=L+5
            ENDIF
         ENDIF
         consecutive=.false.
      ELSE
         consecutive=.true.
         IF ( I .eq. N )  THEN
            write(TERM,'(''-'',I4)') INDX(I)
            KEY=KEY(1:L)//TERM
            L=L+5
         ENDIF
      ENDIF
   ENDDO
!
!     REMOVING BLANK SPACES
!
   J=0
   DO I=1,L
      IF ( KEY(I:I) .ne. ' ') THEN
         J=J+1
         TMPKEY(J:J)=KEY(I:I)
      ENDIF
   ENDDO
   L=J
   KEY=''
   KEY(1:L)=TMPKEY(1:L)

   RETURN
END SUBROUTINE COLSPEC




