beta=1
f=0.05
for propTyped in 0.5 0.2 0.1 0.05 0.02 0.01 ; do #0.2 0.1 0.05 0.01
	#i=$(awk -F, -v af=$f '$3 == af && $1 > 1200 && $1 < 1500 {print $1}' test_2_alleleFreq.csv | head -n 1)
	
	echo $i $f 
	
	#python ../argwas/ARGWAS.py --task associate --out /data/ARGWAS/experiments/freq${f}_index${i}_beta${beta}_withNoise --tree_file /data/ARGWAS/experiments/test_2.trees --name one_variant --ass_method both --pty_sim_method fixed --pty_fixed_variant_indeces $i --pty_fixed_betas $beta --pty_sd_envNoise 1 &
	
	python ../argwas/ARGWAS.py --task associate --out /data/ARGWAS/experiments/freq${f}_indexUntyped_beta${beta}_propTyped${propTyped}_withNoise --tree_file /data/ARGWAS/experiments/test_2.trees --name one_variant_untyped --ass_method both --pty_sim_method singleUntyped --single_variant_af 0.4 --prop_typed_variants $propTyped --pty_fixed_betas $beta --pty_sd_envNoise 1 &

	python ../argwas/ARGWAS.py --task associate --out /data/ARGWAS/experiments/freq${f}_indexUntyped_beta${beta}_propTyped${propTyped} --tree_file /data/ARGWAS/experiments/test_2.trees --name one_variant_untyped --ass_method both --pty_sim_method singleUntyped --single_variant_af 0.4 --prop_typed_variants $propTyped --pty_fixed_betas $beta --pty_sd_envNoise 0 &
	
done




