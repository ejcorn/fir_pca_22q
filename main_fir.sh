#!/bin/bash
set -euxo pipefail

[ -z $BASH ] || shopt -s expand_aliases
alias BEGINCOMMENT="if [ ]; then"
alias ENDCOMMENT="fi"

#################
### Load data ###
#################

# get rid of all of these variables
ZDIM=1
XCP='xcp_6p_noFilter'
#XCP='xcp_36p_despike'
ATLAS=Schaefer
NPARC=200
SCAN='ID'
LAB=''

ROOT='CPCA_'${SCAN}${ATLAS}${NPARC}'Z'${ZDIM}${XCP}${LAB}
LOGDIR='/cbica/home/cornblae/logfiles/'
BASEDIR='/cbica/home/cornblae/ecornblath/fir_pca_22q/'
MATPATH=/cbica/software/external/matlab/R2017A/bin/matlab
MATPATHFAST=/cbica/software/external/matlab/R2014B/bin/matlab # here put an earlier version of matlab if possible ... makes it easier to get jobs
PYPATH=python
RPATH=/usr/lib64/R/bin/Rscript
QUEUE=all.q
LONGJOB=h_rt=99:00:00,s_rt=99:00:00
MASTERDIR=$BASEDIR'results/'$ROOT

mkdir -p $MASTERDIR		# recursively create general results folder and output folder

cd $BASEDIR'jobs'		# change to directory containing all shell scripts

####################
### Process data ###
####################

#qsub -N "PROCESSDEMO" -o ${LOGDIR} -e ${LOGDIR} -v BD=$BASEDIR,RP=$RPATH ProcessDemo.sh
BEGINCOMMENT
qsub -N "PROCESSIMAGINGDATA" -o ${LOGDIR} -e ${LOGDIR} -v ZDIM=$ZDIM,ATLAS='Schaefer',NPARC=200,SCAN=$SCAN,XCP=$XCP,LAB=$LAB,BD=$BASEDIR,MP=$MATPATH,RP=$RPATH -hold_jid "PROCESSDEMO" ProcessData22q.sh
qsub -N "PROCESSIMAGINGDATA" -o ${LOGDIR} -e ${LOGDIR} -v ZDIM=$ZDIM,ATLAS='HarvardOxford',NPARC=112,SCAN=$SCAN,XCP=$XCP,LAB=$LAB,BD=$BASEDIR,MP=$MATPATH,RP=$RPATH -hold_jid "PROCESSDEMO" ProcessData22q.sh

###########################
### FIR design matrices ###
###########################

qsub -N "PROCESSBEHAVIOR"  -o ${LOGDIR} -e ${LOGDIR} -l h_vmem=1G,s_vmem=1G -q $QUEUE -v D=$ROOT,BD=$BASEDIR,RP=$RPATH -hold_jid "PROCESSIMAGINGDATA" ProcessBehavior.sh

qsub -N "FIR_DESIGN" -o ${LOGDIR} -e ${LOGDIR} -l h_vmem=10.5G,s_vmem=10G -q $QUEUE -v D=$ROOT,BD=$BASEDIR,MP=$MATPATHFAST -hold_jid "PROCESSBEHAVIOR" fir_design.sh

###############################################
### FIR CPCA model and LME group comparison ###
###############################################
ENDCOMMENT
COMP_DESIGN='ThreatNonthreatAllStimuliStratified' # use a regressor that gets threat, non-threat, correct, incorrect
BEGINCOMMENT
# fit CPCA

qsub -N "FIR_CPCA" -o ${LOGDIR} -e ${LOGDIR} -q $QUEUE -v D=$ROOT,BD=$BASEDIR,CD=$COMP_DESIGN,MP=$MATPATHFAST -hold_jid "FIR_DESIGN" fir_cpca.sh

# fit linear mixed effects models to capture CPC responses
qsub -N "LME_FIT" -o ${LOGDIR} -e ${LOGDIR} -q $QUEUE -l h_vmem=2G,s_vmem=2G -v D=$ROOT,BD=$BASEDIR,CD=$COMP_DESIGN,RP=$RPATH -hold_jid "FIR_CPCA" lme_fit.sh
ENDCOMMENT
qsub -N "LME_PLOT" -o ${LOGDIR} -e ${LOGDIR} -q $QUEUE -l h_vmem=2G,s_vmem=2G -v D=$ROOT,BD=$BASEDIR,CD=$COMP_DESIGN,RP=$RPATH -hold_jid "LME_FIT" lme_plot.sh
qsub -N "FIR_COGNITION" -o ${LOGDIR} -e ${LOGDIR} -q $QUEUE -l h_vmem=2G,s_vmem=2G -v D=$ROOT,BD=$BASEDIR,CD=$COMP_DESIGN,RP=$RPATH -hold_jid "FIR_CPCA" fir_pca_cognition.sh
BEGINCOMMENT

# bootstrap CPCA coefficients

qsub -N "FIR_BOOTPREP" -o ${LOGDIR} -e ${LOGDIR} -l h_vmem=10.5G,s_vmem=10G -q $QUEUE -v D=$ROOT,BD=$BASEDIR,MP=$MATPATHFAST,CD=$COMP_DESIGN -hold_jid "FIR_DESIGN" fir_bootprep.sh

for REP in {1..100}; do
	qsub -N "FIR_BOOTRUN_AllSubjects_${REP}" -o ${LOGDIR} -e ${LOGDIR} -t 1-100 -l h_vmem=4.5G,s_vmem=4.5G,tmpfree=4.5G -q $QUEUE -v D=$ROOT,BD=$BASEDIR,PY=$PYPATH,CD=$COMP_DESIGN,GRP="AllSubjects",REP=$REP -hold_jid "FIR_BOOTPREP" fir_boot.sh		
done

qsub -N "FIR_BOOTANALYZE" -o ${LOGDIR} -e ${LOGDIR} -l h_vmem=10.5G,s_vmem=10G -q $QUEUE -v D=$ROOT,BD=$BASEDIR,MP=$MATPATH,CD=$COMP_DESIGN -hold_jid "FIR_BOOTRUN_AllSubjects_*","FIR_CPCA" fir_bootanalyze.sh

# test whether PNC and 22q coefficients separately will fit the other sample vs. group coefficients, in bootstrapped samples
ENDCOMMENT
qsub -N "FIR_COEFF_COMPARE" -o ${LOGDIR} -e ${LOGDIR} -l h_vmem=10.5G,s_vmem=10G -q $QUEUE -v D=$ROOT,BD=$BASEDIR,MP=$MATPATH,CD=$COMP_DESIGN,RP=$RPATH -hold_jid "FIR_CPCA" fir_coeff_compare.sh

###################################################################################
### Compare PC maps to structural difference maps and white matter connectivity ###
###################################################################################

# see if t1 maps from Sun et al. 2018 align better with components than spun versions of maps with preserved spatial covariance
#qsub -N "MAKE_SPINS" -o ${LOGDIR} -e ${LOGDIR} -q $QUEUE -v D=$ROOT,BD=$BASEDIR,MP=$MATPATHFAST -hold_jid "PROCESSIMAGINGDATA" make_spins.sh
qsub -N "APPLY_SPINS" -o ${LOGDIR} -e ${LOGDIR} -q $QUEUE -v D=$ROOT,BD=$BASEDIR,MP=$MATPATHFAST -hold_jid "MAKE_SPINS" apply_spins.sh
qsub -N "SPIN_TEST" -o ${LOGDIR} -e ${LOGDIR} -q $QUEUE -v D=$ROOT,BD=$BASEDIR,RP=$RPATH -hold_jid "APPLY_SPINS" spin_test.sh

# white matter connectivity from this dataset
qsub -N "SPIN_TEST" -o ${LOGDIR} -e ${LOGDIR} -q $QUEUE -v D=$ROOT,BD=$BASEDIR,RP=$RPATH -hold_jid "APPLY_SPINS" spin_test.sh
