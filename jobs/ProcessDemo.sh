#!/bin/bash
set -euxo pipefail
cd ${BD}code/process
$RP ProcessDemo.R $BD
