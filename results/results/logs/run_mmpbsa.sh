#!/usr/bin/env bash
set -eo pipefail
source /opt/miniforge/bin/activate mmpbsa
export PATH="/opt/miniforge/envs/mmpbsa/bin:$PATH"
export OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
BASE=/root/roots; RESULTS=$BASE/results/mmpbsa; mkdir -p $RESULTS
DT_PS=100; LIG=UNL; NP=4; CORES=16
echo "=== $(date) cores=$(nproc) NP=$NP ==="

prep() {
  local sys=$1 OUTDIR=$RESULTS/$1
  mkdir -p $OUTDIR
  local TPR=$(find $BASE/harvest/$sys -name '*.tpr' | head -1)
  local XTC=$(find $BASE/harvest/$sys -name '*.xtc' | sort | tail -1)
  local FIRST=$(gmx check -f $XTC 2>&1 | grep -oP 'time\s+\K[0-9.]+' | head -1)
  local START=$(python3 -c "print(int(float('$FIRST')+1000))")
  gmx make_ndx -f $TPR -o $OUTDIR/complex.ndx <<MK >$OUTDIR/make_ndx.log 2>&1
r ${LIG}
"Protein" | r ${LIG}
q
MK
  local GL=$(gmx make_ndx -f $TPR -n $OUTDIR/complex.ndx <<<'q' 2>&1)
  local PIDX=$(echo "$GL" | awk '$2=="Protein"&&$3==":"{print $1;exit}')
  local LIDX=$(echo "$GL" | awk '($2=="UNL"||$2=="r_UNL")&&$3==":"{print $1;exit}')
  echo "$PIDX $LIDX" > $OUTDIR/group_indices.txt
  echo -e "Protein\nProtein_${LIG}" | gmx trjconv -f $XTC -s $TPR -n $OUTDIR/complex.ndx -o $OUTDIR/complex_traj.xtc -b $START -dt $DT_PS -pbc mol -center >$OUTDIR/trjconv.log 2>&1
  cat > $OUTDIR/mmpbsa_gb.in <<'G'
&general
sys_name="ROOTS", startframe=1, endframe=9999999, interval=1, forcefields="charmm", PBRadii=3, temperature=310,
/
&gb
igb=8, saltcon=0.150, ifqnt=0,
/
&decomp
idecomp=1, dec_verbose=1, print_res="within 6"
/
G
  cat > $OUTDIR/mmpbsa_pb.in <<'P'
&general
sys_name="ROOTS", startframe=1, endframe=9999999, interval=1, forcefields="charmm", PBRadii=3, temperature=310,
/
&pb
istrng=0.150, fillratio=4.0, inp=1, radiopt=0, indi=1.0, exdi=80.0,
/
&decomp
idecomp=1, dec_verbose=1, print_res="within 6"
/
P
  echo "[$sys] prepped P=$PIDX L=$LIDX from ${START}ps"
}

runone() {
  local sys=$1 cs=$2 ce=$3
  local OUTDIR=$RESULTS/$sys TPR=$(find $BASE/harvest/$sys -name '*.tpr' | head -1)
  local TOP=$BASE/topo/$sys/topol.top
  read PIDX LIDX < $OUTDIR/group_indices.txt
  cd $OUTDIR
  echo "[$sys] GB cores $cs-$ce ($(date +%H:%M:%S))"
  taskset -c ${cs}-${ce} mpirun -np $NP --bind-to none gmx_MMPBSA MPI -O     -i mmpbsa_gb.in -cs $TPR -ci complex.ndx -cg "$PIDX" "$LIDX"     -ct complex_traj.xtc -cp $TOP -o FINAL_RESULTS_MMPBSA_GB.dat -eo FINAL_RESULTS_MMPBSA_GB.csv     -do FINAL_DECOMP_MMPBSA_GB.dat -deo FINAL_DECOMP_MMPBSA_GB.csv -nogui >gb.log 2>&1
  echo "[$sys] PB ($(date +%H:%M:%S))"
  taskset -c ${cs}-${ce} mpirun -np $NP --bind-to none gmx_MMPBSA MPI -O     -i mmpbsa_pb.in -cs $TPR -ci complex.ndx -cg "$PIDX" "$LIDX"     -ct complex_traj.xtc -cp $TOP -o FINAL_RESULTS_MMPBSA_PB.dat -eo FINAL_RESULTS_MMPBSA_PB.csv     -do FINAL_DECOMP_MMPBSA_PB.dat -deo FINAL_DECOMP_MMPBSA_PB.csv -nogui >pb.log 2>&1
  echo "[$sys] COMPLETE ($(date +%H:%M:%S))"
}

for s in Sys1 Sys2 Sys3 Sys4; do prep $s; done

CPS=$((CORES/4))
runone Sys1 0 $((CPS-1)) &
runone Sys2 $CPS $((2*CPS-1)) &
runone Sys3 $((2*CPS)) $((3*CPS-1)) &
runone Sys4 $((3*CPS)) $((4*CPS-1)) &
wait
echo "=== ALL COMPLETE === $(date)"
