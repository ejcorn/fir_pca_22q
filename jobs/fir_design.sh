#!/bin/bash
set -euxo pipefail

$MP -nodisplay -r "name_root = char('$D'); basedir = char('$BD'); cd(basedir); addpath(genpath('code')); run([basedir,'code/fir/fir_design.m']); exit"
