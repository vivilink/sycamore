import os

os.chdir("/home1/linkv/ARGWAS/hawaiian/plink_files_analysis_chromosomes/chr5/")

pca = open("/home1/linkv/ARGWAS/hawaiian/plink_files_wholeGenome_original/PC100/nh.megaandgda.geno05.hwe1e-6.ld50-5-.8.maf01.rmvChr5.100pcs.smartpca.evec", "r")
pca_new = open("nh.megaandgda.geno05.hwe1e-6.ld50-5-.8.maf01.rmvChr5.100pcs.smartpca_reformatted.evec", "w")

#os.chdir("/home1/linkv/ARGWAS/hawaiian/plink_files_analysis_chromosomes/chr16/")

#pca = open("nh.megaandgda.geno05.hwe1e-6.ld50-5-.8.maf01.rmvChr16.smartpca.evec", "r")
#pca_new = open("nh.megaandgda.geno05.hwe1e-6.ld50-5-.8.maf01.rmvChr16.smartpca_reformatted.evec", "w")

eigvals = pca.readline()

#pca_new.write(header)

for line in pca:
	line = line.rstrip("\n")
	line = ' '.join(line.split())
	ind = line.split(" ")[0].split(":")[0]
	newline = ind + " " + ind 
	for i in line.split(" ")[1:101]:
		newline += " " + i 
	newline += "\n"
	pca_new.write(newline)
	

pca.close()
pca_new.close()
