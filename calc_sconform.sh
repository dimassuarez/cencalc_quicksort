#!/bin/bash
#
# This script helps estimate the 
# conformational entropy from AMBER topology file and trajectory files
#

##########################################################################3
# Programs being used


#CENCALC programs
export LD_LIBRARY_PATH="/opt/dislin/":$LD_LIBRARY_PATH
CENCALC_PATH="/opt/ccmla/"
CENCALC="$CENCALC_PATH/cencalc_ccmla"

# Make sure that get_tor.py has execution permissions
GETTOR="$CENCALC_PATH/get_tor.py -noMet -puck" # Choose torsion selection options here 

# AMBER trajectory analysis software
CPPTRAJ="$AMBERHOME/bin/cpptraj"

# Octave is used for plotting data
OCTAVE=$(which octave | grep -v alias)
if [ -z "$OCTAVE" ]; then echo 'OCTAVE is not available, but recommended!'; DO_PLOT=0; else DO_PLOT=1; fi

#  NCDUMP converts netCDF file to text format. NETCDF library is required (usually bundled with AMBER).
NCDUMP=$(which  ncdump| grep -v alias)
if [ -z "$NCDUMP" ]; then echo 'ncudmp seems to be not available, but needed!';  DO_EXIT=1; fi

#  GNU Parallel builds and runs commands in parallel. 
PARALLEL=$(which  parallel | grep -v alias)
if [ -z "$PARALLEL" ]; then echo 'parallel seems to be not available, but needed!';  DO_EXIT=1; fi

##########################################################################3
# The following environmental variables specify the solute mask, 
# topology file, trajectory data, details of CC-MLA, ...
echo "calc_sconform.sh is running...."
DO_EXIT=0

# CCMLA options
if [ -z "$DO_CC_MLA" ]
then
   DO_CC_MLA="0"
fi
if [ -z "$DO_S2" ]
then
   DO_S2="0"
fi
if [ $DO_CC_MLA -eq 1 ]; then echo "Processing trajectory and doing CC-MLA calcs "; fi
if [ $DO_CC_MLA -eq 2 ]; then echo "Doing CC-MLA calcs (previously processsed trajectory) "; fi

if [ $DO_S2 -eq 1 ]; then echo "Processing trajectory and doing S2 calcs "; fi
if [ $DO_S2 -eq 2 ]; then echo "Doing S2 calcs (previously processsed trajectory) "; fi

# Composite CC-MLA ?
if [ -z "$DO_COMP_CC_MLA" ]
then
   DO_COMP_CC_MLA="0"
fi
if [ $DO_CC_MLA -gt 0 ] && [ $DO_COMP_CC_MLA -eq  1 ]; then echo "Trying composite CC-MLA calcs "; fi

# Cutoff values for CCMLA   (CUTOFF=-1 ---> No cutoff)
if [ -z "$CUTOFF" ]
then
   export CUTOFF="6 8 10"
   echo "Using CUTOFF=$CUTOFF values for CCMLA calcs (if requested)"
else
   echo "Using CUTOFF=$CUTOFF values for CCMLA calcs as predefined"
fi

if [ $DO_CC_MLA -le 1 ] && [ $DO_S2 -le 1 ] 
then

# MASK DEFINITION  
if [ -z "$MOL_MASK" ]
then
   echo "MOL_MASK not defined... then keeping only solute atoms"
   MOL_MASK="!:WAT,Na+,Cl-"
fi
echo "Using MOL_MASK=$MOL_MASK"

# TOPOLOGY FILE (Full path name)
if [ -z "$REFTOP" ]
then
   echo "Topology file not defined, but needed." 
   DO_EXIT="1"
else 
   if [ ! -e  "$REFTOP" ]
   then 
       echo "Topology file $REFTOP does not exist"
       DO_EXIT="1"
   fi
   echo "Using REFTOP=$REFTOP"
fi

# Directory containing trajectory files
if [ -z "$TRAJDIR" ]
then
   echo "Trajectory directory not defined, but needed." 
   DO_EXIT=1
else 
   if [ ! -e  "$TRAJDIR" ]
   then 
       echo "Traj directory $TRAJDIR does not exist"
       if [ $DO_CC_MLA -le 1 ];then  DO_EXIT=1; fi 
   else
      echo "Using TRAJDIR=$TRAJDIR"
   fi
fi
if [ -z "$MDCRD_PREFIX" ]
then
   echo "Prefix of mdcrd files not defined, but needed." 
   DO_EXIT=1
else 
   echo "Assuming MDCRD_PREFIX=$MDCRD_PREFIX"
fi
if [ -z "$MDCRD_SUFFIX" ]
then
   echo "Suffix of mdcrd files not defined, but recommended." 
   MDCRD_SUFFIX=".mdcrd"
   echo "Attempting MDCRD_SUFFIX=$MDCRD_SUFFIX"
else 
   echo "Assuming MDCRD_SUFFIX=$MDCRD_SUFFIX"
fi

fi     # endif of DO_CC_MLA -lt 1

# TEMP Directories
if [  -z "$SCRATCH" ] 
then
   echo "SCRATCH directory not defined, but needed." 
   DO_EXIT=1
else 
   if [ ! -e  "$SCRATCH" ]
   then 
       echo "SCRATCH=$SCRATCH does not exist"
       DO_EXIT=1
   else 
       TT=$(date +%N)
       TMPDIR=${SCRATCH}/TMPDIR_${TT}
       echo "Created $TMPDIR temporal directory"
       echo "Warning: /dev/shm (RAM) will be also used!"
   fi
fi

# Number of processors
if [ -z "$NPROCS" ]
then
   export NPROCS=$(cat /proc/cpuinfo | grep -c processor)
   echo "Using all $NPROCS available processors"
else
   echo "Using $NPROCS processors as predefined "
fi
export OMP_NUM_THREADS=$NPROCS 
export OMP_STACKSIZE="2G"    
 
# The following variables remain normally unchanged

# OFFSET   
if [ -z "$OFFSET" ]
then
   echo "Processing all snapshots"
   OFFSET="1"
else
   echo "Processing snapshots with OFFSET=$OFFSET"
fi

# NINTERVAL
if [ -z "$NINTERVAL" ]
then
   NINTERVAL="20"
   echo "Convergence plot using ${NINTERVAL} points"
else
   echo "Convergence plot using ${NINTERVAL} points"
fi

# Concentration Parameter
if [ -z "$KPARAM" ]
then
   KPARAM="0.50"
fi
echo "Concentration parameter K=$KPARAM"

# Maximum number of conformations (3 recommended)
if [ -z "$MAXNUMCONF" ]
then
   MAXNUMCONF="3"
fi
echo "Maximum number of conformations per torsion=$MAXNUMCONF"

# VERBOSE option for CCMLA
if [ -z "$VERBOSE" ]
then
   VERBOSE=""
else
   VERBOSE="-verbose"
fi

# Cutoffs for Composite Calcs
if [ -z "$CUTS2F1" ]
then
    CUTS2F1="0.01"
fi
if [ -z "$CUTS2F2" ]
then
    CUTS2F2="0.10"
fi

# DISTANCE MATRIX CALCULATION (ON/OFF)
DISTMAT="ON" 

if [ $DO_EXIT -eq 1 ]; then  echo 'Sorry, cannot continue'; exit; fi

##########################################################################3
# Run chicken run 

#WORK directory
WORKDIR=$PWD

if [ "${DO_CC_MLA}" -le 1 ] && [ "${DO_S2}" -le 1 ]  
then

# TRAJECTORY 
i=0
declare -a TRAJ=""
for midterm in $(cd $TRAJDIR; ls ${MDCRD_PREFIX}*${MDCRD_SUFFIX} | sed "s/${MDCRD_PREFIX}//" | sed "s/${MDCRD_SUFFIX}//" | sort -n; cd $WORKDIR)
do
   let "i=$i+1"
   TRAJ["$i"]=$TRAJDIR/${MDCRD_PREFIX}${midterm}${MDCRD_SUFFIX}
done
NTRAJ="$i"

echo "After doing ls $TRAJDIR/${MDCRD_PREFIX}*${MDCRD_SUFFIX} .... we get"
echo "number of trajectory files = $NTRAJ "  

if [ $NTRAJ -eq 0 ]
then 
    echo "No trajectory files"; exit;
fi

# Working in the local scratch may improve efficiency 
mkdir $TMPDIR
cd    $TMPDIR 

# Number of frames in the trajectory file
for ((index=1; index <= NTRAJ; index++))
do

NFRAMES[$index]=$($NCDUMP -h  ${TRAJ[$index]} | grep 'frame = ' | sed -e 's/.*[^0-9]\([0-9]\+\)[^0-9]*$/\1/')

INIT[$index]="1"

echo "Trajectory: ${TRAJ[$index]} Init= ${INIT[$index]} nframes= ${NFRAMES[$index]}  offset= ${OFFSET}" 

done

# Preparing cpptraj input 

$GETTOR -sel ${MOL_MASK/:/} $REFTOP  > temp.in 

for ((index=1; index <= NTRAJ; index++))
do
   echo  "trajin ${TRAJ[$index]}  ${INIT[$index]}  ${NFRAMES[$index]}  10"   >>  matrix.in 
done

sed  -n '/matrix/,$ p' temp.in   >> matrix.in 

for ((index=1; index <= NTRAJ; index++))
do
  echo  "trajin ${TRAJ[$index]}  ${INIT[$index]}  ${NFRAMES[$index]}  ${OFFSET} "    >>  torsion.in 
done

grep -v trajin temp.in | sed '/matrix/,$ d' >> torsion.in
rm -f  temp.in

# Running CPPTRAJ
if [ $NPROCS -gt "1" ]
then 

# Splitting the torsion.in file 
grep -v 'trajin' torsion.in | grep -v '#' > temp.in
nlines=$(cat temp.in | wc -l)
if [ $nlines -gt $OMP_NUM_THREADS ]
then
   let " nlines_per_task = ( $nlines / $OMP_NUM_THREADS ) "
else
   nlines_per_task=1
fi
split -l $nlines_per_task  -d  temp.in   temp_X
rm -f input.torsion 
itor=0
for file in $(ls temp_X*) 
do
  let "itor=$itor+1"
  for ((index=1; index <= NTRAJ; index++))
  do
    echo  "trajin ${TRAJ[$index]}  ${INIT[$index]}  ${NFRAMES[$index]}  ${OFFSET} "    >>  torsion_${itor}.in
  done
  cat $file >> torsion_${itor}.in
  echo " ${CPPTRAJ}  $REFTOP  <  torsion_${itor}.in > torsion_${itor}.out " >>input.torsion 
  rm -f $file 
done 

cat input.torsion | $PARALLEL -t -j $NPROCS
    
else

echo "Running cpptraj torsion/pucker ..."
$CPPTRAJ  $REFTOP < torsion.in > torsion_cpptraj.out 

fi

# CPPTRAJ CPU TIME 
CPU_TOR=$(grep 'Total execution time:' torsion*out | awk '{print SUM+=$5}' | tail -1)

# Processing dihedral data 
ls d????.dat > LISTA

rm -f TASK.sh
echo "Running cencalc_prep individual tasks (time consuming)"
for file in $(cat LISTA)
do
  dihed=${file/.dat/}
  echo  "$PREP -k $KPARAM  -maxconf $MAXNUMCONF -nocut -ag yes -plot yes -wrbigmat ${dihed}.bm ${file} > ${dihed}.out" >> TASK.sh
done

chmod 755 TASK.sh
export OMP_NUM_THREADS=1
cat TASK.sh | $PARALLEL -j $NPROCS

CPU_PREP=$(grep ' CPU-TIME ' d????.out | awk '{print SUM+=$4}' | tail -1)
paste d????.bm | sed 's/\t//g' > BIGMAT.dat
grep -h 'Minima' d????.out  > MINIMA.dat

export OMP_NUM_THREADS=$NPROCS

if [ $DISTMAT == "ON" ]
then

  echo "Running cpptraj matrix job ..."
  $CPPTRAJ   $REFTOP < matrix.in > matrix_cpptraj.out 
  CPU_DISTMAT=$(grep 'Total execution time:' matrix_cpptraj.out | awk '{print $5}' )

  #Preparing CENCALC data
  echo "Running cencalc_prep global ..."
  $PREP -k $KPARAM  -maxconf $MAXNUMCONF -ag yes -plot yes -rdbigmat BIGMAT.dat  d????.dat > cencalc_prep.out 

else

#Preparing CENCALC data
  echo "Running cencalc_prep global ..."
  $PREP -k $KPARAM  -maxconf $MAXNUMCONF -ag yes -plot yes -nocut  -rdbigmat BIGMAT.dat  d????.dat > cencalc_prep.out 

fi

# Copy important stuff in WORKDIR
cp cencalc_prep.out MATRIX.dat MINIMA.dat reduced_dist_matrix.dat distance_matrix.dat \
   matrix_cpptraj.out atoms_in_tor.info torsion.in d*.png    $WORKDIR/

echo "CPU TIME INFO : (s)" >  $WORKDIR/CPU_TIME.dat
echo "===================" >> $WORKDIR/CPU_TIME.dat
echo "Cpptraj-dihedral : $CPU_TOR " >> $WORKDIR/CPU_TIME.dat
echo "Cpptraj-distmat  : $CPU_DISTMAT " >> $WORKDIR/CPU_TIME.dat
echo "Cencalc-Prep     : $CPU_PREP " >> $WORKDIR/CPU_TIME.dat

#Running cencalc 
echo "Running cencalc S1 ..."

NFRAMES=$(wc -l MATRIX.dat | awk '{print $1}')
let "NSTEP = ${NFRAMES} / ${NINTERVAL}"
$CENCALC -s1  -c -1 -data MATRIX.dat  -ns $NSTEP $NFRAMES $NSTEP -t s1_plot.tab > s1.out

CPU_S1=$(grep ' Total ' s1.out  | grep CPU | awk '{print $4}') 
echo "Cencalc-S1       : $CPU_S1 " >> $WORKDIR/CPU_TIME.dat

cp  s1_plot.tab  s1.out $WORKDIR/ 

cd $WORKDIR/

rm -r -f $TMPDIR 

fi 

if [ "$DO_CC_MLA" -eq 1 ]  || [ "$DO_CC_MLA" -eq 2 ]
then

   if [ ! -e MATRIX.dat ];  then  echo "DO_CC_MLA=1,2 but MATRIX.dat does not exist"; exit; fi
   NFRAMES=$(wc -l MATRIX.dat | awk '{print $1}')
   let "NSTEP = ${NFRAMES} / ${NINTERVAL}"

   TT=$(date +%N)
   TMPSHM=/dev/shm/TMPDIR_${TT}
   mkdir $TMPSHM
   cp MATRIX.dat $TMPSHM/

   rm -f ccmla_plot.tab
   touch ccmla_plot.tab 
   cuttoff_plot=""

   if [ ! -e s1_plot.tab ]; then $CENCALC -s1 -c -1 -data $TMPSHM/MATRIX.dat  -ns $NSTEP $NFRAMES $NSTEP -t s1_plot.tab > s1.out; fi

   if [ $DO_COMP_CC_MLA -eq  0 ]
   then 

        for cutoff in $(echo $CUTOFF)
        do
            echo "Running cencalc CCMLA with cutoff=${cutoff} "
            $CENCALC $VERBOSE -ccmla  -c ${cutoff} -data $TMPSHM/MATRIX.dat   -dist  reduced_dist_matrix.dat \
             -ns $NSTEP  $NFRAMES  $NSTEP  -t s_ccmla_${cutoff}.tab >  s_ccmla_${cutoff}.out
             paste ccmla_plot.tab s_ccmla_${cutoff}.tab  > tmp; mv tmp ccmla_plot.tab
             cutoff_plot="${cutoff_plot} ${cutoff}"
             CPU_CC=$(grep ' Total ' s_ccmla_${cutoff}.out  | grep CPU | awk '{print $4}') 
             echo "Cencalc-CCMLA ${cutoff}  : $CPU_CC " >> CPU_TIME.dat
        done
        sed -i 's/\t/ /g' ccmla_plot.tab

   else

       touch ccmla_plot.tab 
       cuttoff_plot=""
       for cutoff in $(echo $CUTOFF)
       do
           echo "Running cencalc composite CCMLA cutoff=${cutoff} "
           $CENCALC $VERBOSE  -ccmla  -c ${cutoff} -data $TMPSHM/MATRIX.dat   -dist  reduced_dist_matrix.dat \
             -s2filt $CUTS2F1  $CUTS2F2  -prs2mat s_ccmla_comp.s2mat \
            -ns $NSTEP  $NFRAMES  $NSTEP  -t s_ccmla_comp_${cutoff}.tab >  s_ccmla_comp_${cutoff}.out
           paste ccmla_plot.tab s_ccmla_comp_${cutoff}.tab  > tmp; mv tmp ccmla_plot.tab
           cutoff_plot="${cutoff_plot} ${cutoff}"
           CPU_CC=$(grep ' Total ' s_ccmla_comp_${cutoff}.out  | grep CPU | awk '{print $4}') 
           echo "Cencalc-CCMLA comp $CUTS2F1  $CUTS2F2  ${cutoff} : $CPU_CC " >> CPU_TIME.dat
       done
       sed -i 's/\t/ /g' ccmla_plot.tab

fi
     
$OCTAVE -q <<EOF
A=load('s1_plot.tab');
s1=A(:,2);
isnap=A(:,1);
clear A;

cutoff=[ $cutoff_plot  ];

nplot=length(cutoff);

B=load('ccmla_plot.tab');
[n,m]=size(B);

ncorr=(m/nplot);
j=0;
scc=zeros(n,nplot);
for i=[ncorr:ncorr:m]
  j=j+1;
  scc(:,j)=B(:,i);
end

T=0.300;
%Representacion grafica 
clf();
h=figure(1);

plot (isnap,-T*s1,'-o','Linewidth',2,'Markersize',4)
ymin=min(-T*s1);
ymax=max(-T*s1);
legend('s1');
hold
txt=['S_1'];
composite=${DO_COMP_CC_MLA};

for j=[1:nplot]
  plot (isnap,-T*scc(:,j),'-o','Linewidth',2,'Markersize',4)
  if composite == 1
  txt=[txt;['CCMLA Comp. r_c=',num2str(cutoff(j))]];
  else
  txt=[txt;['CCMLA r_c=',num2str(cutoff(j))]];
  end
  ymin_corr=min(-T*scc(:,j));
  ymax_corr=max(-T*scc(:,j));
  ymin=min([ymin,ymin_corr]);
  ymax=max([ymax,ymax_corr]);
end

grid on;
W = 8; H = 6;
set(h,'PaperUnits','inches')
set(h,'PaperOrientation','portrait');
set(h,'PaperSize',[H,W])
set(h,'PaperPosition',[0,0,W,H])
set(gca,'Fontname','Times')
set(gca,'Fontsize',18)
legend(txt);
xlim([0 max(isnap)])
ylim([ymin*1.1  ymax*0.9 ])
xlabel('# Snap ','Fontsize',18)
% \305 for Angstroms 
ylabel(['-TS_{conform} kcal/mol'],'Fontsize',18)
print(h,'s_ccmla.png','-dpng','-color')

if nplot > 1 

if composite == 1 
   fid=fopen('s_ccmla_composite.tab','w');
else
   fid=fopen('s_ccmla.tab','w');
end
fmtA=' %s ';
fmtB=' %i ';
for j=[1:nplot]
   if  composite == 1
   fmtA=[fmtA,' CCMLA Comp r_c=%6.2f '];
   else
   fmtA=[fmtA,' CCMLA r_c=%6.2f '];
   end
   fmtB=[fmtB,' %10.4f '];
end
fmtA=[fmtA,'\n'];
fmtB=[fmtB,'\n'];
fprintf(fid,fmtA,'# NumSnap',cutoff);
for i=[1:n]
   fprintf(fid,fmtB,isnap(i),scc(i,:))
end
fclose(fid);

end

EOF

rm -r -f $TMPSHM 

fi


if [ "$DO_S2" -eq 1 ]  || [ "$DO_S2" -eq 2 ]
then

   if [ ! -e MATRIX.dat ];  then  echo "DO_S2=1,2 but MATRIX.dat does not exist"; exit; fi
   NFRAMES=$(wc -l MATRIX.dat | awk '{print $1}')
   let "NSTEP = ${NFRAMES} / ${NINTERVAL}"

   TT=$(date +%N)
   TMPSHM=/dev/shm/TMPDIR_${TT}
   mkdir $TMPSHM
   cp MATRIX.dat $TMPSHM/

   rm -f s2_plot.tab
   touch s2_plot.tab 
   cuttoff_plot=""

   if [ ! -e s1_plot.tab ]; then $CENCALC -s1 -c -1 -data $TMPSHM/MATRIX.dat  -ns $NSTEP $NFRAMES $NSTEP -t s1_plot.tab > s1.out; fi

        for cutoff in $(echo $CUTOFF)
        do
            echo "Running cencalc S2 with cutoff=${cutoff} "
            $CENCALC $VERBOSE -s2 -shuffle  -c ${cutoff} -data $TMPSHM/MATRIX.dat   -dist  reduced_dist_matrix.dat \
             -ns $NSTEP  $NFRAMES  $NSTEP  -t s2_${cutoff}.tab >  s2_${cutoff}.out
             paste s2_plot.tab s2_${cutoff}.tab  > tmp; mv tmp s2_plot.tab
             cutoff_plot="${cutoff_plot} ${cutoff}"
             CPU_S2=$(grep ' Total ' s2_${cutoff}.out  | grep CPU | awk '{print $4}') 
             echo "Cencalc-S2 ${cutoff}  : $CPU_S2 " >> CPU_TIME.dat
        done
        sed -i 's/\t/ /g' s2_plot.tab
     
$OCTAVE -q <<EOF
A=load('s1_plot.tab');
s1=A(:,2);
isnap=A(:,1);
clear A;

cutoff=[ $cutoff_plot  ];

nplot=length(cutoff);

B=load('s2_plot.tab');
[n,m]=size(B);

ncorr=(m/nplot);
j=0;
scc=zeros(n,nplot);
for i=[ncorr:ncorr:m]
  j=j+1;
  scc(:,j)=B(:,i);
end

T=0.300;
%Representacion grafica 
clf();
h=figure(1);

plot (isnap,-T*s1,'-o','Linewidth',2,'Markersize',4)
ymin=min(-T*s1);
ymax=max(-T*s1);
legend('s1');
hold
txt=['S_1'];

for j=[1:nplot]
  plot (isnap,-T*scc(:,j),'-o','Linewidth',2,'Markersize',4)
  txt=[txt;['S2 r_c=',num2str(cutoff(j))]];
  ymin_corr=min(-T*scc(:,j));
  ymax_corr=max(-T*scc(:,j));
  ymin=min([ymin,ymin_corr]);
  ymax=max([ymax,ymax_corr]);
end

grid on;
W = 8; H = 6;
set(h,'PaperUnits','inches')
set(h,'PaperOrientation','portrait');
set(h,'PaperSize',[H,W])
set(h,'PaperPosition',[0,0,W,H])
set(gca,'Fontname','Times')
set(gca,'Fontsize',18)
legend(txt);
xlim([0 max(isnap)])
ylim([ymin*1.1  ymax*0.9 ])
xlabel('# Snap ','Fontsize',18)
% \305 for Angstroms 
ylabel(['-TS_{conform} kcal/mol'],'Fontsize',18)
print(h,'s2.png','-dpng','-color')

if nplot > 1 

fid=fopen('s_s2.tab','w');
fmtA=' %s ';
fmtB=' %i ';
for j=[1:nplot]
   fmtA=[fmtA,' S2 r_c=%6.2f '];
   fmtB=[fmtB,' %10.4f '];
end
fmtA=[fmtA,'\n'];
fmtB=[fmtB,'\n'];
fprintf(fid,fmtA,'# NumSnap',cutoff);
for i=[1:n]
   fprintf(fid,fmtB,isnap(i),scc(i,:))
end
fclose(fid);

end

EOF




rm -r -f $TMPSHM 

fi

#=======================================================================================
 

