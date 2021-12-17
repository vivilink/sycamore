~/git/relate/bin/Relate --haps ${1}_variants.haps --sample ${1}_inds.sample --mode All --output ${1} --mutation_rate 1.25e-8 --effectiveN 500 --map genetic_map_GRCh37_chr1.map

~/git/relate/bin/RelateFileFormats --input ${1} --output ${1}  --mode ConvertToTreeSequence
