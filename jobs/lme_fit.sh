#!/bin/bash
set -euxo pipefail
cd ${BD}code/fir
$RP fir_pca_lme_stim.R $D $BD $CD
