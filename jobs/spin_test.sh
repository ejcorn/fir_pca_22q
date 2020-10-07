#!/bin/bash
set -euxo pipefail
cd ${BD}code/t1
$RP plot_t1_spins.R $D $BD $CD
