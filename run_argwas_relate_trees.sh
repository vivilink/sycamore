beta=1
af=0.4
for propTyped in 0.7 0.5 0.2 0.1 0.05 0.02 0.01 
do 		
	python ../argwas/ARGWAS.py --task associate --out /data/ARGWAS/experiments/freq${af}_indexUntyped_beta${beta}_propTyped${propTyped}_withNoise --tree_file /data/ARGWAS/experiments/relate_test_2_propTyped${propTyped}/test_2_propTyped${propTyped}.trees --tree_file_simulated /data/ARGWAS/experiments/test_2.trees --variants_file /data/ARGWAS/experiments/relate_test_2_propTyped${propTyped}/test_2_propTyped${propTyped}_filtered_sample_variants.csv  --name one_variant_untyped --ass_method both --pty_sim_method singleUntyped --single_variant_af $af --pty_fixed_betas $beta  --pty_sd_envNoise 1 &

	python ../argwas/ARGWAS.py --task associate --out /data/ARGWAS/experiments/freq${af}_indexUntyped_beta${beta}_propTyped${propTyped} --tree_file /data/ARGWAS/experiments/relate_test_2_propTyped${propTyped}/test_2_propTyped${propTyped}.trees --tree_file_simulated /data/ARGWAS/experiments/test_2.trees --variants_file /data/ARGWAS/experiments/relate_test_2_propTyped${propTyped}/test_2_propTyped${propTyped}_filtered_sample_variants.csv --name one_variant_untyped --ass_method both --pty_sim_method singleUntyped --single_variant_af $af --pty_fixed_betas $beta  --pty_sd_envNoise 0 &
done

#need to simulate typed variant if propTyped=1
#propTyped=1

#python ../argwas/ARGWAS.py --task associate --out /data/ARGWAS/experiments/freq${af}_indexUntyped_beta${beta}_propTyped${propTyped} --tree_file /data/ARGWAS/experiments/test_2.trees --name one_variant_untyped --ass_method both --pty_sim_method singleUntyped --single_variant_af $af --pty_fixed_betas $beta  --pty_sd_envNoise 0 --min_allele_freq 0.01


# simulated tree
#propTyped=simulated

#python ../argwas/ARGWAS.py --task associate --out /data/ARGWAS/experiments/freq${af}_indexUntyped_beta${beta}_propTyped${propTyped} --tree_file /data/ARGWAS/experiments/test_2.trees --name one_variant_untyped --ass_method both --pty_sim_method singleTyped --single_variant_af $af --pty_fixed_betas $beta  --pty_sd_envNoise 0 --min_allele_freq 0.01

#--variants_file /data/ARGWAS/experiments/test_2_propTyped${propTyped}_filtered_sample_variants.csv
#--variants_file /data/ARGWAS/experiments/test_2_propTyped${propTyped}_filtered_sample_variants.csv


