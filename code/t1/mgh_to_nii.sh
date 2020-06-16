#!/bin/bash
set -euxo pipefail
cd data/Sunetal2018_Fig1Maps

MGH=( $(find *.mgh -type f))
for i in $(seq ${#MGH[@]}); do	
	mri_convert -it mgh -ot nii ${MGH[i]} ${MGH[i]%.mgh}.nii
done