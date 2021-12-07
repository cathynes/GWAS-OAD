### scripot for GWAS
for chr in {1..22}; do
        
echo "#!/bin/sh 

/home/users/tcheandj/phewas/plink2 \
  --covar $PATH_to_covar_files/$covar_files.txt  \
  --covar-name PC1-PC10 Age_At_Recruitment Sex BSA  (covariate names) \
  --glm hide-covar cols=+a1freq,+machr2  \
  --mac 10 \
  --mach-r2-filter 0.4 2.0 \
  --memory 42000 \
  --pheno-name $phenotype \
  --out $path_to_output/result.chr${chr} \
  --pfile /oak/stanford/projects/ukbb/genotypes/pgen_app13721_v3/ukb_imp_chr${chr}_v3.mac1 \
  --remove /oak/stanford/groups/jpriest/catherine/data/UKBB_subject.to.be.exclude.txt \
  --pheno /$PATH_to_phenotype_files/$phenotype.txt \
  --threads 6 " > $path_to_script/script/gwas.$phenotype.chr${chr}.sh 


chmod u=rwx $path_to_script/script/gwas.$phenotype.chr${chr}.sh 

sbatch -p normal,owners,jpriest --qos=normal -J job_names.chr${chr} -t 24:00:00 -N 1 -n 1 \
-o $PATH_to_log_file/job_error.chr${chr}.log \
-e $PATH_to_log_file/job_log.chr${chr}.err --mem=42000 \
--open-mode=append $path_to_script/script/gwas.$phenotype.chr${chr}.sh ;

done

 ## save the script MI_GWAS_UKBB.sh 

### running the script
./MI_GWAS_UKBB.sh