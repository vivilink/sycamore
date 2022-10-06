#!/bin/bash

gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1  --HEreg --grm ${1} --pheno ${1}_phenotypes.phen --out ${1}_GRM_covariance_tests --threads 2 > ${1}_tmp.out

#extract correct lines, replace white space with tab and remove possibly introduced double tabs
sed -n '2,4p' ${1}_GRM_covariance_tests.HEreg | unexpand -a | tr -s '\t' > ${1}_HE-CP_result.txt
sed -n '7,9p' ${1}_GRM_covariance_tests.HEreg | unexpand -a | tr -s '\t' > ${1}_HE-SD_result.txt



