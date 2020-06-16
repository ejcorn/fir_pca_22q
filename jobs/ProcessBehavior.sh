#!/bin/bash
set -euxo pipefail

cd $BD'code/behavior/'
$RP extract_TRlocked_responses_v2.R $D $BD