beta=1
af=0.4
for propTyped in 0.7 0.5 0.2 0.1 0.05 0.02 0.01 
do 		

	dir=/home/vivian/postdoc_USC/AIM/experiments
	treeFile=$dir/relate_test_2_propTyped${propTyped}/test_2_propTyped${propTyped}.trees
	treeFileSim=$dir/test_2.trees
	variantsFile=$dir/relate_test_2_propTyped${propTyped}/test_2_propTyped${propTyped}_filtered_sample_variants.csv
	outPrefix=$dir/freq${af}_indexUntyped_beta${beta}_propTyped${propTyped}_withNoise

	python ~/git/argwas/ARGWAS.py --task associate --out $outPrefix --tree_file $treeFile --tree_file_simulated $treeFileSim --variants_file $variantsFile  --name one_variant_untyped --ass_method both --pty_sim_method allelicHetero --allelic_hetero_file allelicHetero_instructions.txt  --pty_sd_envNoise 1 &
	
	python ~/git/argwas/ARGWAS.py --task associate --out $outPrefix --tree_file $treeFile --tree_file_simulated $treeFileSim --variants_file $variantsFile  --name one_variant_untyped --ass_method both --pty_sim_method allelicHetero --allelic_hetero_file $dir/allelicHetero_instructions.txt  --pty_sd_envNoise 1
	
	outPrefix=$dir/freq${af}_indexUntyped_beta${beta}_propTyped${propTyped}

	python ~/git/argwas/ARGWAS.py --task associate --out $outPrefix --tree_file $treeFile --tree_file_simulated $treeFileSim --variants_file $variantsFile --name one_variant_untyped --ass_method both --pty_sim_method singleUntyped --single_variant_af $af --pty_fixed_betas $beta  --pty_sd_envNoise 0 &
	
	python ~/git/argwas/ARGWAS.py --task associate --out $outPrefix --tree_file $treeFile --tree_file_simulated $treeFileSim --variants_file $variantsFile --name one_variant_untyped --ass_method both --pty_sim_method allelicHetero --allelic_hetero_file $dir/allelicHetero_instructions.txt  --pty_sd_envNoise 0
done

#need to simulate typed variant if propTyped=1
#propTyped=1

#python ../argwas/ARGWAS.py --task associate --out /data/ARGWAS/experiments/freq${af}_indexUntyped_beta${beta}_propTyped${propTyped} --tree_file /data/ARGWAS/experiments/test_2.trees --name one_variant_untyped --ass_method both --pty_sim_method singleUntyped --single_variant_af $af --pty_fixed_betas $beta  --pty_sd_envNoise 0 --min_allele_freq 0.01


# simulated tree
#propTyped=simulated

#python ../argwas/ARGWAS.py --task associate --out /data/ARGWAS/experiments/freq${af}_indexUntyped_beta${beta}_propTyped${propTyped} --tree_file /data/ARGWAS/experiments/test_2.trees --name one_variant_untyped --ass_method both --pty_sim_method singleTyped --single_variant_af $af --pty_fixed_betas $beta  --pty_sd_envNoise 0 --min_allele_freq 0.01

#--variants_file /data/ARGWAS/experiments/test_2_propTyped${propTyped}_filtered_sample_variants.csv
#--variants_file /data/ARGWAS/experiments/test_2_propTyped${propTyped}_filtered_sample_variants.csv


python ~/git/argwas/ARGWAS.py --task associate --out $outPrefix --tree_file $treeFile --tree_file_simulated $treeFileSim --variants_file $variantsFile  --name one_variant_untyped --ass_method both --pty_sim_method allelicHetero --allelic_hetero_file $dir/allelicHetero_instructions.txt  --pty_sd_envNoise 1
