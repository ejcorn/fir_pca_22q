#!/bin/bash
set -euxo pipefail

$MP -nodisplay -r "name_root = char('$D'); basedir = char('$BD'); component_design = char('$CD'); fin=$FIN; st = $ST; cd(basedir); addpath(genpath('code')); run([basedir,'code/fir/fir_pca_compare_coeffs_ci.m']); exit"
$RP code/fir/fir_varexp_specificity.R $D $BD $CD $FIN $ST
