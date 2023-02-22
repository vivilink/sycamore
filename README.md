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


Analysis descriptions
-------------------------

*simulated ARGs for power analysis*

I simulated 300 random ARG. On the cluster they are located here: /home1/linkv/ARGWAS/power_sims/tree_files/stdpopsim/normal_trees. I downsampled the variants with an allele frequency of at least 10% to 20% typed variants and estimated Relate trees from these variants here: /home1/linkv/ARGWAS/power_sims/tree_files/stdpopsim/normal_trees/minFreq0.01. I also estimated Relate trees from all variants here: ~/ARGWAS/power_sims/tree_files/stdpopsim/normal_trees/allVariants_relate

*random phenotypes and assoiation tests for cutoff*

The null simulations consisting of random phenotypes and association tests for the true trees and all variants are located here: ~/ARGWAS/simulations_cutoff/stdpopsim/N2K/diploid/eGRM_GRM/true_trees/window_based. For the Relate trees and downsampled variants they are located here: /home1/linkv/ARGWAS/simulations_cutoff/stdpopsim/N2K/diploid/eGRM_GRM/relate_trees/window_based. For the Relate trees estimated from all variants they are located here: ~/ARGWAS/simulations_cutoff/stdpopsim/N2K/diploid/eGRM_GRM/relate_trees_allVariants/window_based/. The directories are further divided into directories 5k and 10k, which contain the association test results with corresponding window sizes.

The R script used to calculate the cutoff values is: R_scripts/1_calculate_significance_cutoffs.R

*association tests for power analysis*

The main directory is /home1/linkv/ARGWAS/power_sims/stdpopsim/. It is further separated by true trees, relate trees (these are estimated from 20% of the true trees' variants) and relate_trees_allVariants (these are estimated from all variants of the true trees), and phenotypes with a single causal variant (oneVariant) and phenotypes with allelic heterogeneity (oneRegion). For the paper, I'm using the results in the eGRM_GRM and window_based folders. For allelic heterogeneity, the next distinction is the causal window size (10k or 5k) and the testing window size (tested10k or tested5k). For the single variant case, there are tests for testing window size (tested10k or tested5k) and causal variant allele frequency = 0.02 (rareVariant) and frequency = 0.2 (commonVariant). I started running allelic heterogeneity for causal window size 10k and testing window size 5k, but the results are not complete and also not based on the newest way of running things (the relate trees and true trees do not use the same phenotypes although they should)

The R script used to calculate the association power is 2a_calculate_power.R, which calls 2b_calculate_power_one_experiment.R. It uses the output of R_scripts/1_calculate_significance_cutoffs.R. The results can be plotted with R_scripts/3_plot_with_error_bars_aH.R and R_scripts/3_plot_with_error_bars_oneVariant.R.

The folders with association results called "ownPhenotypes" are old results where the phenotypes were simulated separately for the true trees and the relate trees, and the number of causal variants is thus not consistent

*population structure*

The ARGs of two populations that I used in the end is simulated here: ~/ARGWAS/simulations_two_populations/tree_files/no_migration_highNe_splitTime10k. This is also where the global eGRMs are located. The association tests that correct for population structure are located here: ~/ARGWAS/simulations_two_populations/association_tests/eGRM_GRM/true_trees/window_based/5k/with_strat_correction_highNe_splitTime10k, and the association tests that do not correct for populations structure are here: ~/ARGWAS/simulations_two_populations/association_tests/eGRM_GRM/true_trees/window_based/5k/no_strat_correction_highNe_splitTime10k

The R script to plot the results is R_scripts/population_structure.R.

*CREBRF with REML*

I first estimated piecewise eGRMs for every RELATE tree of the genome using script ~/ARGWAS/hawaiian/run_egrm.sh. The global eGRMs are located here: ~/ARGWAS/hawaiian/global_grm/. I then combined the piecewise eGRMs, except for those of chromosome 5, to a global eGRM using script R_scripts/combine_egrms.R. I ran association testing for all parts of chromosome 5 with script ~/ARGWAS/hawaiian/run_association.sh, but not all finished. For the paper, I extracted the region around CREBRF into tree file ~/ARGWAS/hawaiian/chr5.part-04.CREBRF.trees, and tested this region for association using ~/ARGWAS/hawaiian/run_association_CREBRF.sh. 

*CREBRF with GWAS*

The original plink files with genotypes around the CREBRF gene and 10 principle components calculated based on the whole genome, which were shared with me, are here:/home1/linkv/ARGWAS/hawaiian/plink_files_copy/ (this is a copy for safety). The analysis for the paper are here: ~/ARGWAS/hawaiian/plink_files/. I first added the standardized phenotypes that I calculated in the REML analysis to the .fam file with other_scripts/add_phenotype_to_fam.py. I then ran /home1/linkv/ARGWAS/hawaiian/plink_files/run_plink.sh to transform the .fam and .bed files into .ped and .map files, and to perform the GWAS association tests.


Other
-----------

- In order to run several of the tasks, you need to provide two ARG files, one with parameter "tree_file" and the other with "tree_file_simulated". in addition to variant files (which are produced when you simulate phenotypes). The tree_file is the main ARG that will be tested for association testing. In the case where you have a RELATE tree, the tree_file is the RELATE tree. The "tree_file_simulated" becomes necessary in the case where you simulate phenotypes based on non-typed variants. Since the RELATE trees only contain typed variants, there needs to be a way to access the origial, non-typed variants. Tree_file_simulated should always be accompanied by the "variant_file" that is produced by task"downsampleVariantsWriteShapeit". If you are working with the true trees and all variants, you can provide the same simulated tree file with both ARG paramters, and the variant file that is produced by task "simulate".
- For the paper, I only use covariance types "eGRM" and "GRM". The "scaled" covariance was the first one implemented and might not have all functionalities. It is calculated based on the TMRCA and the assumption of brownian motion.
- For the paper, I only use the mixed model association results produced with REML. However, the program can also run the Haseman-Elston algorithm of GCTA. HE is faster but has lower power than REML at small sample sizes, e.g. 1000 diploids.
- In the code, the names "argwas" and "AIM" are used interchangeably for our new method that uses local eGRMs to test for associations
- Currently, you need to provide the absolute path for all files (or the relative path from the software folder), please feel free to change that
