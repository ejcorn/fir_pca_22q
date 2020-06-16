#!/bin/bash
set -euxo pipefail

cd $BD
$PY code/fir/boot_fir_pca.py $BD $D $CD $GRP $REP

#python code/fir/boot_fir_pca.py '/cbica/home/cornblae/ecornblath/brain_states_22q/' 'ScanIDSchaefer200Z1_22q_PreQC_XCP36pdespike_us_100reps' 'ThreatNonthreatAllStimuliStratified' 'PNC' 1
