#!/bin/bash

./egrm/bin/trees2egrm --output ${2} ${3} --c-extension --output-format gcta ${1} 

./egrm/bin/trees2egrm --output ${2} ${3} --c-extension --output-format numpy ${1} 

#~/git/argwas/egrm/bin/trees2egrm --output "/data/ARGWAS/experiments_N500/nullSimulation/REML/haploid/test" --c-extension  --skip-first-tree --output-format gcta --verbose /data/ARGWAS/experiments_N500/nullSimulation/REML/haploid/null_propTyped0.1_envNoise1_null_ploidy1_focal.trees
