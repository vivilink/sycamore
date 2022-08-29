#!/bin/bash

gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --reml --grm ${1} --pheno ${1}_phenotypes.phen --out ${1}_REML --threads 2 --reml-maxit 500 > ${1}_tmp.out

#extract correct lines, replace white space with tab and remove possibly introduced double tabs
grep "Pval" ${1}_REML.hsq | cut -f 2 > ${1}_REML_result.txt



#~/git/argwas/gcta_1.93.3beta2/gcta64 --reml --grm freq0.4_indexUntyped_beta1_propTyped0.7_aH_GRM_covariance --pheno freq0.4_indexUntyped_beta1_propTyped0.7_aH_phenotypes.phen --out REML_test --threads 2
