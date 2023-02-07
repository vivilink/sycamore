import os

os.chdir("/data/ARGWAS/hawaiians/CREBRF/plink_GWAS")
inds_pheno = {}

# get transformed BMI

transformed_phenos = open("chr5.part-00_eGRM_phenotypes.phen", "r")
for line in transformed_phenos:
	line = line.strip('\n')
	ind = line.split(" ")[1]
	pheno = line.split(" ")[2]
	inds_pheno[ind] = [pheno]

transformed_phenos.close()


# get sex

orig_phenos = open("crebrf.region.n5383.MEGA+GDA_converted.csv", "r")
for line in orig_phenos:
	line = line.strip('\n')
	ind = line.split(" ")[0]
	sex = line.split(" ")[1]
	entry = inds_pheno.get(ind)
	if entry:
		inds_pheno[ind].append(sex)
	else:
		print("individual " + ind + " does not have a phenotype")
	
orig_phenos.close()	


# add to fam
	
fam = open("crebrf.region.n5383.MEGA+GDA.fam", "r")
fam_new = open("crebrf.region.n5383.MEGA+GDA_withPheno.fam", "w")
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

	
	
