#!/bin/bash

./egrm/bin/trees2egrm --output ${2} ${3} --c-extension --skip-first-tree --output-format gcta ${1} --verbose

./egrm/bin/trees2egrm --output ${2} ${3} --c-extension --skip-first-tree --output-format numpy ${1} --verbose

#~/git/argwas/egrm/bin/trees2egrm --output "/data/ARGWAS/experiments_N500/nullSimulation/REML/haploid/test" --c-extension  --skip-first-tree --output-format gcta --verbose /data/ARGWAS/experiments_N500/nullSimulation/REML/haploid/null_propTyped0.1_envNoise1_null_ploidy1_focal.trees
