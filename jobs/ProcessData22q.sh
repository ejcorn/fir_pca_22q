#!/bin/bash
set -euxo pipefail

$MP -nodisplay -r "addpath(genpath(fullfile('$BD','code'))); basedir = char('$BD'); cd(basedir); zdim = [$ZDIM]; scan = char('$SCAN'); atlasName = char('$ATLAS'); atlasScale = [$NPARC]; extralabel = char('$LAB'); XCP_folder = char('$XCP'); run([basedir,'code/process/ProcessData22q.m']); exit" 

D='CPCA_'${SCAN}${ATLAS}${NPARC}'Z'${ZDIM}${XCP}
cd ${BD}/code/process
$RP samplecharacteristics.R $BD $D