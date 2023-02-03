argwas
=======

Creating conda environment
-------------------------

*Build environment from file*

The list of packages installed in my conda environment used to produce all results can be found here in this file: ARGWAS-package-list.txt. It should thus be possible to recreate my environment with the following command:

    conda create -n myenv --file ARGWAS-package-list.txt
    
*Build environment from scratch*

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

Typical command to run associations with simulated phenotypes and allelic heterogeneity:

    python ARGWAS.py --task simulatePhenotypes --out power_simulation --tree_file relate_ARG.trees --tree_file_simulated simulated_ARG.trees --variants_file simulated_ARG_variants.csv --ass_method GWAS AIM:eGRM AIM:GRM --AIM_method HE REML --pty_sim_method oneRegion --pty_prop_causal_mutations 0.1 --causal_region_coordinates 49500000 49505000 --pty_h_squared 0.02 --pty_sd_beta_causal_mutations standardized --ploidy 2 --skip_first_tree --min_allele_freq 0.01 --seed 1 --ass_window_size 5000 --trees_interval_start 49000000 --simulate_phenotypes
    
Run association with simulated phenotypes and two populations, correcting for stratification with PCA:

    python ARGWAS.py --task associate --out example --tree_file simulated_ARG.trees --tree_file_simulated simulated_ARG.trees  --ass_method AIM:eGRM --AIM_method REML --pty_sim_method null --pty_sd_envNoise 1 --ploidy 2 --seed 1 --min_allele_freq 0 --simulate_phenotypes --population_structure prefix_global_eGRM --population_structure_pca_num_eigenvectors 20 --add_1_to_half_of_inds


*simulate*

This task simulates ARGs using either msprime or stdpopsim. Always simulate haploids, the haplotypes are assigned to individuals in the association task.

Simulate one population with stdpopsim:

    python ARGWAS.py --task simulate --out example --ploidy 1 --sequence_length 30000000.0 --sim_tree_simulator stdPopsim --N 1000 --trees_interval 49000000 50000000 --seed 1

Simulate two populations with msprime:

    python ARGWAS.py --task simulate --mu 1e-08 --out example --ploidy 1 --recomb_rate 1e-08 --sequence_length 100000.0 --sim_tree_simulator msprime --sim_two_populations --N 2000 --seed 1 --split_time 10000 


*simulateMoreMutations*

This task reads in an ARG and adds more mutations in a specified region of the ARG. It then outputs a new ARG.

*downsampleVariantsWriteShapeit*

This task reads a tree, downsamples the variants to mimic a genotyping array and then outputs the variant information in the correct format to be read by Relate to infer ARGs.

    python ARGWAS.py --task downsampleVariantsWriteShapeit --out downsampled --tree_file simulated_ARG.trees --tree_file_simulated simulated_ARG.trees --min_allele_freq 0.01 --ploidy 2 --prop_typed_variants 0.2 --seed 1

    ./relate/bin/Relate --haps downsampled_variants.haps --sample downsampled_inds.sample --mode All --output downsampled_relate --mutation_rate 1.25e-8 --effectiveN 2000 --map genetic_map_GRCh37_chr1.map

    ./relate/bin/RelateFileFormats --input downsampled_relate --output downsampled_relate  --mode ConvertToTreeSequence



*ARGStatistics*

This task takes an ARG and outputs a file with information about the ARG, such as the marginal tree coordinates, number of variants per tree, ...

    python ARGWAS.py --task ARGStatistics --out example --tree_file example.trees --tree_file_simulated example.trees


*simulatePhenotypes*

This task takes an ARG and simulates phenotypes without running any association test. If the same seed is used, the phenotypes simulated under task 'associate' and 'simulatePhenotypes' should be identical.

    python ARGWAS.py --task simulatePhenotypes --out example --tree_file relate_tree.trees --tree_file_simulated simulated_tree.trees --variants_file simulated_tree_variants.csv --ass_method GWAS AIM:eGRM AIM:GRM --AIM_method HE REML --pty_sim_method oneRegion --pty_prop_causal_mutations 0.02 --causal_region_coordinates 49500000 49505000 --pty_h_squared 0.1 --pty_sd_beta_causal_mutations standardized --ploidy 2 --skip_first_tree --min_allele_freq 0.01 --seed 1 --ass_window_size 5000 --trees_interval_start 49000000 --simulate_phenotypes

*impute*

This task takes two ARGs, one for the population of interest and one for a reference panel. It then runs impute2 to infer the missing genotypes in the population of interest using the reference panel.

*getTreeAtPosition*

This task takes an ARG and extracts the marginal tree overlapping a genomic position, and then writes it to file.

*covarianceCorrelations*

This task takes pickle files containing covariance matrices that can be calculated and written to pickle with task 'associate'. It calculates the correlation between these covariance matrices.
