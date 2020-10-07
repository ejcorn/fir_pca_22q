#!/bin/bash
set -euxo pipefail

$MP -nodisplay -r "name_root = char('$D'); basedir = char('$BD'); cd(basedir); addpath(genpath('code')); run([basedir,'code/t1/t1_spin_pcs.m']); exit"
