import os

os.chdir("/data/ARGWAS/hawaiians/CREBRF/plink_GWAS")

pca = open("All_MEGA_GWAS.and.CIDR_GDA_GWAS.nh.intersect1M.drop5k.chr1-22.asVCF.rmvoutoforder.geno05.hwe1e-6.ld50-5-.8.prune.in.maf01.smartpca.evec", "r")
pca_new = open("All_MEGA_GWAS.and.CIDR_GDA_GWAS.nh.intersect1M.drop5k.chr1-22.asVCF.rmvoutoforder.geno05.hwe1e-6.ld50-5-.8.prune.in.maf01.smartpca_reformatted.evec", "w")

eigvals = pca.readline()

#pca_new.write(header)

for line in pca:
	line = line.rstrip("\n")
	line = ' '.join(line.split())
	ind = line.split(" ")[0].split(":")[0]
	newline = ind + " " + ind 
	for i in line.split(" ")[1:11]:
		newline += " " + i 
	newline += "\n"
	pca_new.write(newline)
	

pca.close()
pca_new.close()
