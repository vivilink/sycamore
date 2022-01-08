#!/bin/bash

Rscript create_gcta_GRM.R ${1}

gcta_1.93.3beta2/gcta64 --HEreg --grm ${1}_GRM_covariance --pheno ${1}_phenotypes.phen --out ${1}_GRM_covariance_tests --threads 2 > ${1}_tmp.out

#extract correct lines, replace white space with tab and remove possibly introduced double tabs
sed -n '2,4p' ${1}_GRM_covariance_tests.HEreg | unexpand -a | tr -s '\t' > ${1}_HE-CP_result.txt
sed -n '7,9p' ${1}_GRM_covariance_tests.HEreg | unexpand -a | tr -s '\t' > ${1}_HE-SD_result.txt



