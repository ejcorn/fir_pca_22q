#!/bin/bash
set -euxo pipefail
cd ${BD}code/fir
$RP plot_fir_pca_lme_stim.R $D $BD $CD
$RP fir_pca_lme_coef_tab.R $D $BD $CD