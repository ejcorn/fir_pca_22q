#!/bin/bash
set -euxo pipefail

rsync -avzh /data/jag/bassett-lab/xiaosong/22q/qsirecon/*/*/*.mat cornblae@cbica-cluster.uphs.upenn.edu:/cbica/home/cornblae/ecornblath/brain_states_22q/data/qsiprep/