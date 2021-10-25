#!/bin/bash

Rscript /data/ARGWAS/argwas/create_gcta_GRM.R ${1}

/data/ARGWAS/gcta_1.93.2beta/gcta64 --HEreg --grm /data/ARGWAS/argwas/GRM_covariance_${1} --pheno /data/ARGWAS/argwas/phenotypes_${1}.phen --out /data/ARGWAS/argwas/GRM_covariance_tests_${1} --threads 8 > out_${1}.txt

#extract correct lines, replace white space with tab and remove possibly introduced double tabs
sed -n '2,4p' GRM_covariance_tests_${1}.HEreg | unexpand -a | tr -s '\t' > HE-CP_${1}_result.txt
sed -n '7,9p' GRM_covariance_tests_${1}.HEreg | unexpand -a | tr -s '\t' > HE-SD_${1}_result.txt

