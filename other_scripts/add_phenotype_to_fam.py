import os

os.chdir("/home1/linkv/ARGWAS/hawaiian/plink_files_analysis_chromosomes")
inds_pheno = {}

# get transformed BMI

transformed_phenos = open("~/ARGWAS/hawaiian/all_chr5_for_review/chr5.part-04_chunk28_pca20_phenotypes.phen", "r")
for line in transformed_phenos:
	line = line.strip('\n')
	ind = line.split(" ")[1]
	pheno = line.split(" ")[2]
	inds_pheno[ind] = [pheno]

transformed_phenos.close()


# get sex

orig_phenos = open("~/ARGWAS/hawaiian/MEGA_GDA5_bmi_age_sex.050922.csv", "r")
for line in orig_phenos:
    line = line.strip('\n')
    if len(line.split(",")) == 4:
        ind = line.split(",")[0]
        sex = line.split(",")[1]
        entry = inds_pheno.get(ind)
        if entry:
            inds_pheno[ind].append(sex)
        else:
            print("individual " + ind + " does not have a phenotype")
    else:
        print("line does not have 4 entries", line)
	
orig_phenos.close()	


# add to fam
	
fam = open("~/ARGWAS/hawaiian/plink_files_wholeGenome_original/All_MEGA_GWAS.and.CIDR_GDA_GWAS.nh.intersect1M.drop5k.chr1-22.asVCF.rmvoutoforder.geno05.hwe1e-6.ld50-5-.8.prune.in.maf01.fam", "r")
fam_new = open("whole_genome_withPheno.fam", "w")
for line in fam:
	line = line.strip('\n')
	ind = line.split(" ")[0]
	first = line.split(" ")[0:4]
	newline =""
	for i in first:
		newline += i + " "
	entry = inds_pheno.get(ind)
	if entry:
		# with sex
		#newline += entry[1] + " " + entry[0] + "\n"
		#no sex
		newline += "0 " + entry[0] + "\n"
	else:
		newline += "0 -9\n"
	fam_new.write(newline)

fam_new.close()

	
	
