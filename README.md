sycamore
=======

This repository contains the code used for the paper entitled [Tree-based QTL mapping with expected local genetic relatedness matrices (Link et al. 2023)](https://www.biorxiv.org/content/10.1101/2023.04.07.536093v1). It identifies quantitative trait loci by testing the eGRM, a measure of local relatedness based on the ancestral recombination graph, for association with the phenotypes using a variance-components framework. 

The python code simulates or reads in an ancestral recombination graph and phenotype data (using [stdpopsim](https://github.com/popsim-consortium/stdpopsim), [msprime](https://tskit.dev/msprime/docs/stable/intro.html) and [tskit](https://tskit.dev/)), calculates local eGRMs (using [egrm](https://github.com/vivilink/egrm) based on [Fan et al. 2022](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9118131/)), and tests the local eGRMs for association with the phenotypes (using [GCTA](https://yanglab.westlake.edu.cn/software/gcta/#Overview)). The R scripts plot the results.

Setting up sycamore
-------------------------
* Clone repository

* Clone and install [egrm](https://github.com/vivilink/egrm) with pip in the correct conda environment using the instructions in the egrm readme

* Download [GCTA](https://yanglab.westlake.edu.cn/software/gcta/#Overview) into sycamore repository https://yanglab.westlake.edu.cn/software/gcta/#Download

Creating conda environment
-------------------------

*Build environment from file*

The list of packages installed in my conda environment used to produce all results can be found here in this file: sycamore-package-list.txt. It should thus be possible to recreate my environment with the following command:

    conda env create -n env_sycamore --file sycamore-package-list.txt
    
*Build environment from scratch*

Build environment

    conda env create -f argwas_environment.yml

Install glimix-core

    conda install -c conda-forge glimix-core

Install tskit

    conda install -c conda-forge tskit

<!---
Install limix_lmm

    conda install -c conda-forge limix
-->

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

Install stdpopsim

    conda install -c conda-forge stdpopsim


Install pandas_plink

    conda install -c conda-forge pandas-plink or pip install pandas-plink


Install plinkFile library in R if you want to write covariance matrices in .grm.bin format for GCTA

    install.packages("plinkFile")

Parameters
-----------------

Run python ./ARGWAS.py --help to see a list of all parameters with descriptions

Tasks
-------------------------

*associate*

This task runs all association tests, including GWAS, local REML with eGRM and local REML with GRM. You can either simulate phenotypes or read in existing phenotype files. 

Typical command to run associations with simulated phenotypes and allelic heterogeneity:

    python ARGWAS.py --task simulatePhenotypes --out $dir/power_simulation --tree_file $dir/relate_ARG.trees --tree_file_simulated $dir/simulated_ARG.trees --variants_file $dir/simulated_ARG_variants.csv --ass_method GWAS AIM:eGRM AIM:GRM --AIM_method HE REML --pty_sim_method oneRegion --pty_prop_causal_mutations 0.1 --causal_region_coordinates 49500000 49505000 --pty_h_squared 0.02 --pty_sd_beta_causal_mutations freq_dependant --ploidy 2 --skip_first_tree --min_allele_freq 0.01 --seed 1 --ass_window_size 5000 --trees_interval_start 49000000 --simulate_phenotypes
    
Run association with simulated phenotypes and two populations, correcting for stratification with PCA:

    python ARGWAS.py --task associate --out $dir/example --tree_file $dir/simulated_ARG.trees --tree_file_simulated $dir/simulated_ARG.trees  --ass_method AIM:eGRM --AIM_method REML --pty_sim_method null --pty_sd_envNoise 1 --ploidy 2 --seed 1 --min_allele_freq 0 --simulate_phenotypes --population_structure prefix_global_eGRM --population_structure_pca_num_eigenvectors 20 --add_1_to_half_of_inds


*simulate*

This task simulates ARGs using either msprime or stdpopsim. Always simulate haploids, the haplotypes are assigned to individuals in the association task. If the positions are integers it can happen that more than one variants fall on the same position. In that case, the position of the second variant is set to position + 1.

Simulate one population with stdpopsim:

    python ARGWAS.py --task simulate --out $dir/example --ploidy 1 --sim_tree_simulator stdPopsim --N 1000 --trees_interval 49000000 50000000 --seed 1

Simulate two populations with msprime:

    python ARGWAS.py --task simulate --mu 1e-08 --out $dir/example --ploidy 1 --recomb_rate 1e-08 --sim_tree_simulator msprime --sim_two_populations --N 2000 --seed 1 --split_time 10000 


*simulateMoreMutations*

This task reads in an ARG and adds more mutations in a specified region of the ARG. It then outputs a new ARG.

*downsampleVariantsWriteShapeit*

This task reads a tree, downsamples the variants to mimic a genotyping array and then outputs the variant information in the correct format to be read by Relate to infer ARGs.

    python ARGWAS.py --task downsampleVariantsWriteShapeit --out $dir/downsampled --tree_file $dir/simulated_ARG.trees --tree_file_simulated $dir/simulated_ARG.trees --min_allele_freq 0.01 --ploidy 2 --prop_typed_variants 0.2 --seed 1

    ./relate/bin/Relate --haps downsampled_variants.haps --sample downsampled_inds.sample --mode All --output downsampled_relate --mutation_rate 1.25e-8 --effectiveN 2000 --map genetic_map_GRCh37_chr1.map

    ./relate/bin/RelateFileFormats --input downsampled_relate --output downsampled_relate  --mode ConvertToTreeSequence



*ARGStatistics*

This task takes an ARG and outputs a file with information about the ARG, such as the marginal tree coordinates, number of variants per tree, ...

    python ARGWAS.py --task ARGStatistics --out $dir/example --tree_file $dir/example.trees --tree_file_simulated $dir/example.trees


*simulatePhenotypes*

This task takes an ARG and simulates phenotypes without running any association test. If the same seed is used, the phenotypes simulated under task 'associate' and 'simulatePhenotypes' should be identical.

    python ARGWAS.py --task simulatePhenotypes --out $dir/example --tree_file $dir/relate_tree.trees --tree_file_simulated $dir/simulated_tree.trees --variants_file $dir/simulated_tree_variants.csv --ass_method GWAS AIM:eGRM AIM:GRM --AIM_method HE REML --pty_sim_method oneRegion --pty_prop_causal_mutations 0.02 --causal_region_coordinates 49500000 49505000 --pty_h_squared 0.1 --pty_sd_beta_causal_mutations freq_dependant --ploidy 2 --skip_first_tree --min_allele_freq 0.01 --seed 1 --ass_window_size 5000 --trees_interval_start 49000000 --simulate_phenotypes

*impute*

This task takes two ARGs, one for the population of interest and one for a reference panel. It then runs impute2 to infer the missing genotypes in the population of interest using the reference panel.

*getTreeAtPosition*

This task takes an ARG and extracts the marginal tree overlapping a genomic position, and then writes it to file.

*covarianceCorrelations*

This task takes pickle files containing covariance matrices that can be calculated and written to pickle with task 'associate'. It calculates the correlation between these covariance matrices.


Power Analysis commands
-------------------------

This section contains more details about how to use the software to replicate the power simulations in our paper.

*simulated ARGs for power analysis*

We simulated 300 random ARGs. This produces a tskit object with suffix ".trees" and a file with information about all variants contained within the ARG with suffix "_sample_variants.csv".

    python ./ARGWAS.py --task simulate --N 2000 --out $dir/sim_rep1 --trees_interval 49000000 50000000 --ploidy 1 --seed 1 

We downsampled the variants with an allele frequency of at least 1% to 20% of "typed" variants. This produces a file with information about all variants contained within the ARG with suffix "_sample_variants.csv", similarly to task simulate, but it will contain updated information, e.g. whether a variant is typed or not.

    python ./ARGWAS.py --task downsampleVariantsWriteShapeit --out $dir/sim_rep1_propTyped0.2_minAF0.01 --tree_file $dir/sim_rep1.trees --tree_file_simulated $dir/sim_rep1.trees --min_allele_freq 0.01 --ploidy 2 --prop_typed_variants 0.2 --seed 1

We then estimated Relate trees from all typed variants: 

    ./Relate --haps sim_rep1_propTyped0.2_minAF0.01_variants.haps --sample sim_rep1_propTyped0.2_minAF0.01_inds.sample --mode All --output sim_rep1_propTyped0.2_minAF0.01_relate --mutation_rate 1.25e-8 --effectiveN 2000 --map genetic_map_GRCh37_chr1.map

    ./RelateFileFormats --input sim_rep1_propTyped0.2_minAF0.01_relate --output sim_rep1_propTyped0.2_minAF0.01_relate  --mode ConvertToTreeSequence

*random phenotypes and assoiation tests for cutoff*

We tested the simulated ARGs and all variants for association with random phenotypes in windows of 5k and 10k. We did this for the true simulated trees and the downsampled trees reestimated with Relate. This is the command for window size 5k and the Relate trees:

    python ./ARGWAS.py --task associate --out $dir/sim_rep1_propTyped0.2_minAF0.01 --tree_file $dir/sim_rep1_propTyped0.2_minAF0.01_relate.trees --tree_file_simulated $dir/sim_rep1.trees  --ass_method GWAS AIM:eGRM AIM:GRM --AIM_method HE REML --pty_sim_method null --pty_sd_envNoise 1 --ploidy 2 --seed 1 --skip_first_tree --ass_window_size 5000 --simulate_phenotypes --trees_interval_start 49000000

The R script used to calculate the cutoff values is: R_scripts/1_calculate_significance_cutoffs.R

*association tests for power analysis*

We simulated phenotypes for the individuals of 200 simulated ARGs. We did this for the true trees and the relate trees. In order to still have access to the This is the command we used for alleleic heterogeneity with causal variants in one 5kb window and a local heritability of 0.2: 

    python ./ARGWAS.py --task simulatePhenotypes --out power_sims_rep1 --tree_file $dir/sim_rep1_propTyped0.2_minAF0.01_relate.trees --tree_file_simulated $dir/sim_rep1.trees --variants_file $dir/sim_rep1_propTyped0.2_minAF0.01_sample_variants.csv  --pty_sim_method oneRegion --pty_prop_causal_mutations 0.2 --causal_region_coordinates 49500000 49505000 --pty_h_squared 0.2 --pty_sd_beta_causal_mutations freq_dependant --ploidy 2 --skip_first_tree --min_allele_freq 0.01  --trees_interval_start 49000000 --simulate_phenotypes --seed 1

We then used these phenotypes to test the variants and trees for association. The trees were tested in windows of 5kb or 10kb. Here is the command for 5kb:

    python ./ARGWAS.py --task associate --out power_sims_rep1 --tree_file $dir/sim_rep1_propTyped0.2_minAF0.01_relate.trees --tree_file_simulated $dir/sim_rep1.trees --variants_file $dir/propTyped0.2_minAF0.01_sample_variants.csv  --ass_method GWAS AIM:eGRM AIM:GRM --AIM_method HE REML --ploidy 2 --skip_first_tree --ass_window_size 5000 --trees_interval_start 49000000 --seed 1 --pheno_file $dir/power_sims_rep1.phen

The R script used to calculate the association power is 2a_calculate_power.R, which calls 2b_calculate_power_one_experiment.R. It uses the output of R_scripts/1_calculate_significance_cutoffs.R. The results can be plotted with R_scripts/3_plot_with_error_bars_aH.R and R_scripts/3_plot_with_error_bars_oneVariant.R.

Analyzing BMI data
--------------------

We first create a LOCO (leave-one-chromosome-out) eGRM for each chromosome $chr with output directory $dir:
   
    Rscript sycamore/R_scripts/combine_egrms.R list_trees_LOCO_chr${chr}.txt global_egrm_LOCO_chr${chr} $dir
    awk '{print 0 "\t" $2}' global_egrm_LOCO_chr${chr}.grm.id > global_egrm_LOCO_chr${chr}.grm.id_tmp # replace family ID column with '0' for all (I have had issues with GCTA convergence when all individuals were in separate families)
    mv global_egrm_LOCO_chr${chr}.grm.id_tmp global_egrm_LOCO_chr${chr}.grm.id

We then standardize the BMI data according to McCaw et al. (2019). This entails regressing BMI on sex and age (BMI ~ sex + age + age^2) and removing individuals who's residual is more than 6 standard deviations away from the mean. Outliers are defined separately for each sex with a separate regression. 

We regress out the LOCO eGRM with GCTA for each chromosome $chr:

    gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 --reml --reml-pred-rand --out BLUP_chr${chr} --pheno standardized_phenotypes.phen --grm global_egrm_LOCO_chr${chr} --thread-num 10

We extract the residuals from the GCTA output:

    awk '{print $1,$2,$6}' BLUP_chr${chr}.indi.blp > residuals_chr${chr}.txt

We test the residuals for association with the ARG $treedir/$tree:

    python ./ARGWAS.py --task associate --tree_file $treedir/$tree --out $dir/example --ass_method AIM:eGRM --AIM_method REML --ploidy 2 --min_allele_freq 0 --pheno_file residuals_chr${chr}.txt --relate_sample_names relate.sample --ass_window_size 5000 --skip_first_tree

You need to provide the sample names of the tree leaves with parameter 'relate_sample_names', and these will be matched to the sample names in the phenotype file. Individuals with missing phenotypes, or individuals who are missing from the ARG, are removed. 

In order to parallelize the task, we split each ARG $tree into smaller ARGs of length 1Mb with the following command:

    python ./ARGWAS.py --task makeTreeChunks --tree_file $tree --out $dir/$example --chunk_size 1000000 --skip_first_tree


If you do not want to use the two-step approach, here are the parameters to correct for population structure:

- --population_structure_matrix: provide prefix for global (e.g. LOCO) genetic relatedness matrix in GCTA's binary format
- --population_structure_pca_num_eigenvectors: if population structure matrix is provided, this parameter will not include the matrix directly, but instead use GCTA to calculate the specified number of PCs and include those
- --do_all_stratification_correction: if the population structure matrix is provided and you also provided population_structure_pca_num_eigenvectors, then you can choose to include both the matrix and the PCs in the model with this parameter


Other
-----------

- In order to run several of the tasks, you need to provide two ARG files, one with parameter "tree_file" and the other with "tree_file_simulated". The tree_file is the main ARG that will be tested for association testing. In the case where you have a Relate tree, the tree_file is the Relate tree. The "tree_file_simulated" becomes necessary in the case where you simulate phenotypes based on non-typed variants. Since the Relate trees only contain typed variants, there needs to be a way to access the origial, non-typed variants. Tree_file_simulated should always be accompanied by the "variant_file" that is produced by task "downsampleVariantsWriteShapeit". If you are working with the true trees and all variants, you can provide the same simulated tree file with both ARG paramters, and the variant file that is produced by task "simulate".
- For the paper, we only use covariance types "eGRM" and "GRM". The "scaled" covariance was the first one implemented and might not have all functionalities. It is calculated based on the TMRCA and the assumption of brownian motion.
- For the paper, we only use the mixed model association results produced with REML. However, the program can also run the Haseman-Elston algorithm of GCTA. HE is faster but has lower power than REML at small sample sizes, e.g. 1000 diploids.
- In the code, the names "argwas" and "AIM" are used interchangeably for our new method that calculates local eGRMs and tests them for associations with the phenotype
- Currently, you need to provide the absolute path for all files (or the relative path from the software folder)

Questions
-----------

Please contact linkv@usc.edu with any questions.

Cite us
-----------
Use DOI:10.5281/zenodo.10011427 to cite this software
