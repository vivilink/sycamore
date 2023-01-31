argwas
=======

Creating conda environment
-------------------------

Build environment

    conda env create -f argwas_environment.yml

Install tskit

    conda install -c conda-forge tskit

Install limix_lmm

    conda install -c conda-forge limix

Download and unpack gcta

    https://cnsgenomics.com/software/gcta/#Download, provide correct path with parameter --GCTA

Install statsmodel

    conda install -c conda-forge statsmodels

Install log indenter (I think there is no conda package)

    pip install python_log_indenter

Install tsinfer

    conda install -c conda-forge tsinfer 

Install matplotlib

    conda install -c conda-forge matplotlib
 
Install tsdate

    python3 -m pip install tsdate --user


Install plinkFile library in R if you want to write covariance matrices in .grm.bin format for GCTA

    install.packages("plinkFile")
    
    

Tasks
-------------------------

*associate*

This task runs all association tests, including GWAS, local REML with eGRM and local REML with GRM. You can either simulate phenotypes or read in existing phenotype files. 

*simulate*

This task simulates ARGs using either msprime or stdpopsim

*simulateMoreMutations*

This task reads in an ARG and adds more mutations in a specified region of the ARG. It then outputs a new ARG.

*downsampleVariantsWriteShapeit*

This task reads a tree, downsamples the variants to mimic a genotyping array and then outputs the variant information in the correct format to be read by Relate to infer ARGs.

*ARGStatistics*

This task takes an ARG and outputs a file with information about the ARG, such as the marginal tree coordinates, number of variants per tree, ...

*simulatePhenotypes*

This task takes an ARG and simulates phenotypes without running any association test. If the same seed is used, the phenotypes simulated under task 'associate' and 'simulatePhenotypes' should be identical.

*impute*

This task takes two ARGs, one for the population of interest and one for a reference panel. It then runs impute2 to infer the missing genotypes in the population of interest using the reference panel.

*getTreeAtPosition*

This task takes an ARG and extracts the marginal tree overlapping a genomic position, and then writes it to file.

*covarianceCorrelations*

This task takes pickle files containing covariance matrices that can be calculated and written to pickle with task 'associate'. It calculates the correlation between these covariance matrices.
