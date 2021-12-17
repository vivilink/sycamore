for propTyped in 1 0.7  0.5 0.2 0.1 0.05 0.02 0.01; do #0.5 0.2 0.1 0.05 0.02 0.01

	echo $propTyped
	
	python ../argwas/ARGWAS.py --task downsampleVariants --out /data/ARGWAS/experiments/test_2_propTyped${propTyped} --tree_file /data/ARGWAS/experiments/test_2.trees --prop_typed_variants $propTyped
	
	~/git/relate/bin/Relate --haps test_2_propTyped${propTyped}_variants.haps --sample test_2_propTyped${propTyped}_inds.sample --mode All --output test_2_propTyped${propTyped} --mutation_rate 1.25e-8 --effectiveN 500 --map genetic_map_GRCh37_chr1.map 

	~/git/relate/bin/RelateFileFormats --input test_2_propTyped${propTyped} --output test_2_propTyped${propTyped}  --mode ConvertToTreeSequence 

	#mkdir relate_test_2_propTyped${propTyped}

	mv test_2_propTyped${propTyped}.trees relate_test_2_propTyped${propTyped}
	mv test_2_propTyped${propTyped}.anc relate_test_2_propTyped${propTyped}
	mv test_2_propTyped${propTyped}.mut relate_test_2_propTyped${propTyped}
	mv test_2_propTyped${propTyped}_filtered_sample_variants.csv relate_test_2_propTyped${propTyped}
done
