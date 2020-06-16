#!/bin/bash
set -euxo pipefail

$MP -nodisplay -r "name_root = char('$D'); basedir = char('$BD'); component_design = char('$CD'); cd(basedir); addpath(genpath('code')); run([basedir,'code/fir/fir_pca_compare_coeffs_ci.m']); exit"
$RP code/fir_varexp_specificity.R $D $BD $CD
