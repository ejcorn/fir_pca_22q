#!/bin/bash
set -euxo pipefail

BASEDIR=/cbica/home/cornblae/ecornblath/fir_pca_22q/

echo 'Enter root folder:'
#read ROOT
ROOT=CPCA_IDSchaefer200Z1xcp_36p_despike
#ROOT=CPCA_IDSchaefer200Z1xcp_6p_noFilter

mkdir -p results/$ROOT/
#rsync -avzh cornblae@cbica-cluster.uphs.upenn.edu:/data/tesla-data/ecornblath/brain_states_22q/results/$ROOT/analyses/ results/$ROOT/
#ROOTS=('ScanIDLaus250Z1_22q_100reps' 'ScanIDSchaefer200Z1_22q_100reps' 'ScanIDSchaefer400Z1_22q_100reps' 'ScanIDBrainnetome246Z1_22q_100reps')
#for ROOT in ${ROOTS[@]}; do
rsync -avzh cornblae@cubic-login.uphs.upenn.edu:${BASEDIR}results/$ROOT/analyses/* results/$ROOT/analyses/
#done
#rsync -avzh --exclude "repkmeans*" cornblae@cbica-cluster.uphs.upenn.edu:${BASEDIR}results/$ROOT/clusterAssignments results/$ROOT/

