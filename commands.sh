python ARGWAS.py --task simulate --out test_1 


python ARGWAS.py --task associate --out high_freq_high_beta --tree_file test_2.trees --name one_variant --ass_method GWAS --pty_sim_method fixed --pty_fixed_variant_indeces 1227 --pty_fixed_betas 1
python ARGWAS.py --task associate --out high_freq_high_beta --tree_file test_2.trees --name one_variant --ass_method ARGWAS --pty_sim_method fixed --pty_fixed_variant_indeces 1227 --pty_fixed_betas 1

python ARGWAS.py --task associate --out low_freq_high_beta --tree_file test_2.trees --name one_variant --ass_method GWAS --pty_sim_method fixed --pty_fixed_variant_indeces 1201 --pty_fixed_betas 1
python ARGWAS.py --task associate --out low_freq_high_beta --tree_file test_2.trees --name one_variant --ass_method ARGWAS --pty_sim_method fixed --pty_fixed_variant_indeces 1201 --pty_fixed_betas 1
