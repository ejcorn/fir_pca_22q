#!/bin/bash
set -euxo pipefail

$MP -nodisplay -r "name_root = char('$D'); basedir = char('$BD'); fin=$FIN; st = $ST; component_design = char('$CD'); cd(basedir); addpath(genpath('code')); run([basedir,'code/fir/fir_pca_bootstrap_prep.m']); exit"
