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

Typical command:

    python ARGWAS.py --task simulatePhenotypes --out power_simulation --tree_file relate_ARG.trees --tree_file_simulated simulated_ARG.trees --variants_file simulated_ARG_variants.csv --ass_method GWAS AIM:eGRM AIM:GRM --AIM_method HE REML --pty_sim_method oneRegion --pty_prop_causal_mutations 0.1 --causal_region_coordinates 49500000 49505000 --pty_h_squared 0.02 --pty_sd_beta_causal_mutations standardized --ploidy 2 --skip_first_tree --min_allele_freq 0.01 --seed 1 --ass_window_size 5000 --trees_interval_start 49000000 --simulate_phenotypes

*simulate*

This task simulates ARGs using either msprime or stdpopsim

*simulateMoreMutations*

This task reads in an ARG and adds more mutations in a specified region of the ARG. It then outputs a new ARG.

*downsampleVariantsWriteShapeit*

This task reads a tree, downsamples the variants to mimic a genotyping array and then outputs the variant information in the correct format to be read by Relate to infer ARGs.

    python ARGWAS.py --task downsampleVariantsWriteShapeit --out downsampled --tree_file simulated_ARG.trees --tree_file_simulated simulated_ARG.trees --min_allele_freq 0.01 --ploidy 2 --prop_typed_variants 0.2 --seed 1

    ./relate/bin/Relate --haps downsampled_variants.haps --sample downsampled_inds.sample --mode All --output downsampled_relate --mutation_rate 1.25e-8 --effectiveN 2000 --map genetic_map_GRCh37_chr1.map

    ./relate/bin/RelateFileFormats --input downsampled_relate --output downsampled_relate  --mode ConvertToTreeSequence



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
