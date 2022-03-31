#!/bin/bash

gcta_v1.94.0Beta_linux_kernel_3_x86_64/gcta_v1.94.0Beta_linux_kernel_3_x86_64_static --reml --grm ${1}_GRM_covariance --pheno ${1}_phenotypes.phen --out ${1}_REML --threads 2 > ${1}_tmp.out

#extract correct lines, replace white space with tab and remove possibly introduced double tabs
grep "Pval" ${1}_REML.hsq | cut -f 2 > ${1}_REML_result.txt



#~/git/argwas/gcta_1.93.3beta2/gcta64 --reml --grm freq0.4_indexUntyped_beta1_propTyped0.7_aH_GRM_covariance --pheno freq0.4_indexUntyped_beta1_propTyped0.7_aH_phenotypes.phen --out REML_test --threads 2
