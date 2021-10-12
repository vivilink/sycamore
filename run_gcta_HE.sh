#!/bin/bash

Rscript /data/ARGWAS/argwas/create_gcta_GRM.R

/data/ARGWAS/gcta_1.93.2beta/gcta64 --HEreg --grm /data/ARGWAS/argwas/GRM_covariance --pheno /data/ARGWAS/argwas/phenotypes.phen --out /data/ARGWAS/argwas/GRM_covariance_tests --threads 15 > out.txt

#extract correct lines, replace white space with tab and remove possibly introduced double tabs
sed -n '2,4p' GRM_covariance_tests.HEreg | unexpand -a | tr -s '\t' > HE-CP_result.txt
sed -n '7,9p' GRM_covariance_tests.HEreg | unexpand -a | tr -s '\t' > HE-SD_result.txt

