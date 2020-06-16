#!/bin/bash
set -euxo pipefail

BASEDIR=/cbica/home/cornblae/ecornblath/fir_pca_22q/

echo 'Enter root folder:'
#read ROOT
ROOT=xcp_6p_noFilter
#ROOT=xcp_36p_despike
#ROOT=CPCA_IDSchaefer200Z1xcp_6p_noFilter

rsync -avzh cornblae@cubic-login.uphs.upenn.edu:${BASEDIR}data/sc/*$ROOT* data/sc/
