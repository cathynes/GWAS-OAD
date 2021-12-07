
#!/bin/sh

sbatch -p jpriest,normal,owners --qos=normal -J grm_AOD_inds -t 24:00:00 -N 1 -n 10 -o grm_AOD_inds.log -e grm_AOD_inds.err --mem=62000 \
--wrap="
### computing relatness for the subsset of ukbb data I have using gcta
/scratch/users/tcheandj/gcta_1.93.0beta/gcta64 \
 --bfile /oak/stanford/groups/jpriest/ukb_genotypes/ukb_ALL \
 --keep /oak/stanford/groups/jpriest/catherine/data/all.inds.with.AOD_first_second_set.txt \
 --autosome --maf 0.01 --make-grm --thread-num 10 \
 --out /scratch/users/tcheandj/GCTA_COJO_AOD/grm_AOD_inds"

cd /scratch/users/tcheandj/GCTA_COJO_AOD
 ## get the supset of related I will use a cut off of 0.05 (for identifying 3rd degree relative)
 /scratch/users/tcheandj/gcta_1.93.0beta/gcta64 \
 --grm grm_AOD_inds --grm-singleton 0.05  --out list_inds ## this show that out of the 35059 inds, there are 2109 related inds

## I will now get the list of unrelated inds based on the 0.05 cut off 
 /scratch/users/tcheandj/gcta_1.93.0beta/gcta64 \
 --grm grm_AOD_inds --grm-cutoff 0.05 --make-grm --out unrelateness.list


### estimation of heritability on the entire set of data using GCTA

 ### step1 : extract imputed SNP from the main data

 for chr in {1..22}; do

echo "#!/bin/sh 

#### extract SNP all inds fron the pgen files 
/home/users/tcheandj/phewas/plink2 \
  --pfile /oak/stanford/projects/ukbb/genotypes/pgen_app13721_v3/ukb_imp_chr${chr}_v3.mac1 \
  --keep /scratch/users/tcheandj/subset.ukbb.for.clumping/unrelateness.list.grm.id \
  --mac 10 \
  --mach-r2-filter 0.4 2.0 \
  --memory 42000 \
  --make-bed \
  --out /scratch/users/tcheandj/subset.ukbb.for.clumping/unrelated_white_UKBB_with_AOD.chr${chr} 
" > ../script/extract._chr${chr}.sh
chmod u=rwx ../script/extract._chr${chr}.sh

sbatch -p normal,owners,jpriest --qos=normal -J repli_AOD.chr${chr} -t 24:00:00 -N 1 -n 1 -o log.file/repli_AOD.chr${chr}.log  -e log.file/repli_AOD.chr${chr}.log --mem=42000 --open-mode=append ../script/extract._chr${chr}.sh;
done

 ## ### step2 : calculating segment-based LD score
for chr in {1..22}; do

echo "#!/bin/sh 

## estimate LD for each SNP

 /scratch/users/tcheandj/gcta_1.93.0beta/gcta64 \
 --bfile /scratch/users/tcheandj/subset.ukbb.for.clumping/unrelated_white_UKBB_with_AOD.chr${chr} \
 --ld-score-region 200 \
 --out /scratch/users/tcheandj/subset.ukbb.for.clumping/ldscore_UKBB/ldscore_ukbb.chr${chr} " > ../script/ldscore.chr${chr}.sh;

chmod u=rwx ../script/ldscore.chr${chr}.sh

#sbatch -p normal,owners,jpriest --qos=normal -J ldscore.chr${chr} -t 24:00:00 -N 1 -n 1 -o log.file/ldscore.chr${chr}.log -e log.file/ldscore.chr${chr}.err --mem=42000 --open-mode=append ../script/ldscore.chr${chr}.sh;
done

#### step3 : stratify the SNPs by segment-based LD scores in R
setwd("/scratch/users/tcheandj/subset.ukbb.for.clumping/ldscore_UKBB/")
for (chr in 1:22){
	lds_seg = read.table(paste0("ldscore_ukbb.chr",chr,".score.ld"),header=T,colClasses=c("character",rep("numeric",8)))
	
	### I will remove SNP with MAF<0.0001 as well as indel to make the ldscore estimate easier 
	lds_seg$A1=sapply(strsplit(lds_seg$SNP,"_"),"[",2)
	lds_seg$A2=sapply(strsplit(lds_seg$SNP,"_"),"[",3)
	lds_seg$ncharact=nchar(lds_seg$A1)+nchar(lds_seg$A2)
	### remove SNP with more than 2 character
	lds_seg<- lds_seg[lds_seg$ncharact==2,]

	quartiles=summary(lds_seg$ldscore_SNP)

	lb1 = which(lds_seg$ldscore_SNP <= quartiles[2])
	lb2 = which(lds_seg$ldscore_SNP > quartiles[2] & lds_seg$ldscore_SNP <= quartiles[3])
	lb3 = which(lds_seg$ldscore_SNP > quartiles[3] & lds_seg$ldscore_SNP <= quartiles[5])
	lb4 = which(lds_seg$ldscore_SNP > quartiles[5])

	lb1_snp = lds_seg$SNP[lb1]
	lb2_snp = lds_seg$SNP[lb2]
	lb3_snp = lds_seg$SNP[lb3]
	lb4_snp = lds_seg$SNP[lb4]

	write.table(lb1_snp, paste0("/scratch/users/tcheandj/subset.ukbb.for.clumping/ldscore_UKBB/ldby_segment/snp_group1_chr",chr,".txt"), row.names=F, quote=F, col.names=F)
	write.table(lb2_snp, paste0("/scratch/users/tcheandj/subset.ukbb.for.clumping/ldscore_UKBB/ldby_segment/snp_group2_chr",chr,".txt"), row.names=F, quote=F, col.names=F)
	write.table(lb3_snp, paste0("/scratch/users/tcheandj/subset.ukbb.for.clumping/ldscore_UKBB/ldby_segment/snp_group3_chr",chr,".txt"), row.names=F, quote=F, col.names=F)
	write.table(lb4_snp, paste0("/scratch/users/tcheandj/subset.ukbb.for.clumping/ldscore_UKBB/ldby_segment/snp_group4_chr",chr,".txt"), row.names=F, quote=F, col.names=F)

}

### step 4: making GRMs using SNPs stratified into different groups
#!/bin/sh 
 for chr in {1..22}; do

echo "#!/bin/sh 
#### making GRMs using SNPs stratified into different groups 

for i in 1 2 3 4; do
 /scratch/users/tcheandj/gcta_1.93.0beta/gcta64 \
 --bfile /scratch/users/tcheandj/subset.ukbb.for.clumping/unrelated_white_UKBB_with_AOD.chr${chr} \
 --extract /scratch/users/tcheandj/subset.ukbb.for.clumping/ldscore_UKBB/ldby_segment/snp_group$i_chr$chr.txt \
 --make-grm \
 --out /scratch/users/tcheandj/GCTA_COJO_AOD/GRM_GREML_AOD/grm_snp_group$i_chr$chr.txt;
 done
  " > grm_by_group_AoD.chr${chr}.sh;

chmod u=rwx grm_by_group_AoD.chr${chr}.sh

sbatch -p normal,owners,jpriest --qos=normal -J grm_by_group_AoD.chr${chr} -t 24:00:00 -N 1 -n 1 -o ../log.file/grm_by_group_AoD.chr${chr}.log  -e ../log.file/grm_by_group_AoD.chr${chr}.log --mem=22000 --open-mode=append grm_by_group.chr${chr}.sh;
done

#### step 5: perform REML analysis with multiple GRMs
sbatch -p jpriest --qos=normal -J GREML_AoD -t 36:00:00 -N 1 -n 1 -o GREML_AoD.log -e GREML_AoD.err --mem=52000 \
--wrap="#!/bin/sh

### perform REML analysis with multiple GRMs

 /scratch/users/tcheandj/gcta_1.93.0beta/gcta64 \
 --reml --mgrm /scratch/users/tcheandj/GCTA_COJO_AOD/list.grm_by_group_each_chr.txt \
 --pheno /oak/stanford/groups/jpriest/catherine/data/AAoD_pheno_for_herita_in_GCTA.txt  \
 --thread-num 6 --out /scratch/users/tcheandj/GCTA_COJO_AOD/GREML_AoD_white_ukbb_all_set_combined "



 sbatch -p bigmem --qos=normal -J GREML_AoD -t 24:00:00 -N 1 -n 1 -o GREML_AoD.log -e GREML_AoD.err --mem=1600GB --open-mode=append greml_AoD_ukbb.sh

 ### create the list of GRM by chr 
for chr in {1..22}; do 
 grep  chr${chr}.txt  list.grm_by_group_each_chr.txt > list.grm_by_group_chr${chr}.txt;
done

#### step 5 is also performed by chromosome and the final heritability will be the sum of heritability on each chr


### I will also performed HWE to remove SNP with HWE<1e-05
 
 sbatch -p jpriest --qos=normal -J HWE_estimation -t 36:00:00 -N 1 -n 1 -o HWE_estimation.log -e HWE_estimationD.err --mem=22000 \
--wrap="#!/bin/sh
for chr in {1..22}; do

#### extract SNP all inds fron the pgen files 
/home/users/tcheandj/phewas/plink2 \
--bfile /scratch/users/tcheandj/subset.ukbb.for.clumping/unrelated_white_UKBB_with_AOD.chr${chr} \
--hwe 1e-30 --make-just-bim --out /scratch/users/tcheandj/subset.ukbb.for.clumping/SNP_nonHWEdeviation_subset_AoD_ukbb.chr${chr};
done " 



### next import all the variants with non HWE to used to create group of SNP for grm 
library(data.table)
setwd("/scratch/users/tcheandj/subset.ukbb.for.clumping/ldscore_UKBB/")
for (chr in 1:22) {
	lds_seg = read.table(paste0("ldscore_ukbb.chr",chr,".score.ld"),header=T,colClasses=c("character",rep("numeric",8)))
	
	### I will remove  indel to make the ldscore estimate easier 
	lds_seg$A1=sapply(strsplit(lds_seg$SNP,"_"),"[",2)
	lds_seg$A2=sapply(strsplit(lds_seg$SNP,"_"),"[",3)
	lds_seg$ncharact=nchar(lds_seg$A1)+nchar(lds_seg$A2)

	### remove SNP with more than 2 character
	lds_seg<- lds_seg[lds_seg$ncharact==2,]

	### keep only SNP that do not deviate from HWE
	nonHWE<- fread(paste0("/scratch/users/tcheandj/subset.ukbb.for.clumping/SNP_nonHWEdeviation_subset_AoD_ukbb.chr",chr,".bim"),
						header=F)

	names(nonHWE)<- c("chr","SNP","cm","POS","A1","A2")

	lds_seg <- lds_seg[lds_seg$SNP%in%nonHWE$SNP,]

	lds_seg$MAF<- ifelse(lds_seg$freq<0.5, lds_seg$freq,1-lds_seg$freq)
	summary(lds_seg$MAF)
### cut off on quartile and maf
	quartiles=summary(lds_seg$ldscore_SNP)

	lb1 = which(lds_seg$ldscore_SNP <= quartiles[2])
	## now breakdown by maf
	ld1_maf1 = which(lds_seg$ldscore_SNP <= quartiles[2] & lds_seg$MAF <= 0.0001)
	ld1_maf2 = which(lds_seg$ldscore_SNP <= quartiles[2] & lds_seg$MAF > 0.0001 &  lds_seg$MAF <= 0.001)
	ld1_maf3 = which(lds_seg$ldscore_SNP <= quartiles[2] & lds_seg$MAF > 0.001 &  lds_seg$MAF <= 0.05)
	ld1_maf4 = which(lds_seg$ldscore_SNP <= quartiles[2] & lds_seg$MAF > 0.05 &  lds_seg$MAF <= 0.25)
	ld1_maf5 = which(lds_seg$ldscore_SNP <= quartiles[2] & lds_seg$MAF > 0.25)

	ld2 = which(lds_seg$ldscore_SNP > quartiles[2] & lds_seg$ldscore_SNP < quartiles[3])
	## now breakdown by maf
	ld2_maf1 = which(lds_seg$ldscore_SNP > quartiles[2] & lds_seg$ldscore_SNP <= quartiles[3] & lds_seg$MAF <= 0.0001) 
	ld2_maf2 = which(lds_seg$ldscore_SNP > quartiles[2] & lds_seg$ldscore_SNP <= quartiles[3] & lds_seg$MAF > 0.0001 &  lds_seg$MAF <= 0.001)
	ld2_maf3 = which(lds_seg$ldscore_SNP > quartiles[2] & lds_seg$ldscore_SNP <= quartiles[3] & lds_seg$MAF > 0.001 &  lds_seg$MAF <= 0.05)
	ld2_maf4 = which(lds_seg$ldscore_SNP > quartiles[2] & lds_seg$ldscore_SNP <= quartiles[3] & lds_seg$MAF > 0.05 &  lds_seg$MAF <= 0.35)
	ld2_maf5 = which(lds_seg$ldscore_SNP > quartiles[2] & lds_seg$ldscore_SNP <= quartiles[3] & lds_seg$MAF > 0.35) 

	ld3 = which(lds_seg$ldscore_SNP > quartiles[3] & lds_seg$ldscore_SNP < quartiles[5])
	## now breakdown by maf
	ld3_maf1 = which(lds_seg$ldscore_SNP > quartiles[3] & lds_seg$ldscore_SNP < quartiles[5] & lds_seg$MAF <= 0.0001)
	ld3_maf2 = which(lds_seg$ldscore_SNP > quartiles[3] & lds_seg$ldscore_SNP < quartiles[5] & lds_seg$MAF > 0.0001 &  lds_seg$MAF <= 0.001)
	ld3_maf3 = which(lds_seg$ldscore_SNP > quartiles[3] & lds_seg$ldscore_SNP < quartiles[5] & lds_seg$MAF > 0.001 &  lds_seg$MAF <= 0.05)
	ld3_maf4 = which(lds_seg$ldscore_SNP > quartiles[3] & lds_seg$ldscore_SNP < quartiles[5] & lds_seg$MAF > 0.05 &  lds_seg$MAF <= 0.25)
	ld3_maf5 = which(lds_seg$ldscore_SNP > quartiles[3] & lds_seg$ldscore_SNP < quartiles[5] & lds_seg$MAF > 0.25)
	

	ld4 = which(lds_seg$ldscore_SNP >= quartiles[5])
	## now breakdown by maf
	ld4_maf1 = which(lds_seg$ldscore_SNP >= quartiles[5] & lds_seg$MAF <= 0.0001)
	ld4_maf2 = which(lds_seg$ldscore_SNP >= quartiles[5] & lds_seg$MAF > 0.0001 &  lds_seg$MAF <= 0.001)
	ld4_maf3 = which(lds_seg$ldscore_SNP >= quartiles[5] & lds_seg$MAF > 0.001 &  lds_seg$MAF <= 0.05)
	ld4_maf4 = which(lds_seg$ldscore_SNP >= quartiles[5] & lds_seg$MAF > 0.05 &  lds_seg$MAF <= 0.25)
	ld4_maf5 = which(lds_seg$ldscore_SNP >= quartiles[5] & lds_seg$MAF > 0.25)	

#### within each group I will create different maf group <0.05; 0.05-0.15; 0.16-0.25; 0.26-0.35; >0.35

	#lb1_snp = lds_seg$SNP[lb1]
	lb1_snp_maf1 = lds_seg$SNP[ld1_maf1]
	lb1_snp_maf2 = lds_seg$SNP[ld1_maf2]
	lb1_snp_maf3 = lds_seg$SNP[ld1_maf3]
	lb1_snp_maf4 = lds_seg$SNP[ld1_maf4]
	lb1_snp_maf5 = lds_seg$SNP[ld1_maf5]

	#lb2_snp = lds_seg$SNP[lb2]
	lb2_snp_maf1 = lds_seg$SNP[ld2_maf1]
	lb2_snp_maf2 = lds_seg$SNP[ld2_maf2]
	lb2_snp_maf3 = lds_seg$SNP[ld2_maf3]
	lb2_snp_maf4 = lds_seg$SNP[ld2_maf4]
	lb2_snp_maf5 = lds_seg$SNP[ld2_maf5]

	#lb3_snp = lds_seg$SNP[lb3]
	lb3_snp_maf1 = lds_seg$SNP[ld3_maf1]
	lb3_snp_maf2 = lds_seg$SNP[ld3_maf2]
	lb3_snp_maf3 = lds_seg$SNP[ld3_maf3]
	lb3_snp_maf4 = lds_seg$SNP[ld3_maf4]
	lb3_snp_maf5 = lds_seg$SNP[ld3_maf5]

	#lb4_snp = lds_seg$SNP[lb4]
	lb4_snp_maf1 = lds_seg$SNP[ld4_maf1]
	lb4_snp_maf2 = lds_seg$SNP[ld4_maf2]
	lb4_snp_maf3 = lds_seg$SNP[ld4_maf3]
	lb4_snp_maf4 = lds_seg$SNP[ld4_maf4]
	lb4_snp_maf5 = lds_seg$SNP[ld4_maf5]

	#write.table(lb1_snp, paste0("/scratch/users/tcheandj/subset.ukbb.for.clumping/ldscore_UKBB/ldby_segment/snp_group1_chr",chr,".txt"), row.names=F, quote=F, col.names=F)
	write.table(lb1_snp_maf1, paste0("/scratch/users/tcheandj/subset.ukbb.for.clumping/ldscore_UKBB/ldby_segment/snp_group1_maf1_chr",chr,".txt"), row.names=F, quote=F, col.names=F)
	write.table(lb1_snp_maf2, paste0("/scratch/users/tcheandj/subset.ukbb.for.clumping/ldscore_UKBB/ldby_segment/snp_group1_maf2_chr",chr,".txt"), row.names=F, quote=F, col.names=F)
	write.table(lb1_snp_maf3, paste0("/scratch/users/tcheandj/subset.ukbb.for.clumping/ldscore_UKBB/ldby_segment/snp_group1_maf3_chr",chr,".txt"), row.names=F, quote=F, col.names=F)
	write.table(lb1_snp_maf4, paste0("/scratch/users/tcheandj/subset.ukbb.for.clumping/ldscore_UKBB/ldby_segment/snp_group1_maf4_chr",chr,".txt"), row.names=F, quote=F, col.names=F)
	write.table(lb1_snp_maf5, paste0("/scratch/users/tcheandj/subset.ukbb.for.clumping/ldscore_UKBB/ldby_segment/snp_group1_maf5_chr",chr,".txt"), row.names=F, quote=F, col.names=F)

	#write.table(lb2_snp, paste0("/scratch/users/tcheandj/subset.ukbb.for.clumping/ldscore_UKBB/ldby_segment/snp_group2_chr",chr,".txt"), row.names=F, quote=F, col.names=F)
	write.table(lb2_snp_maf1, paste0("/scratch/users/tcheandj/subset.ukbb.for.clumping/ldscore_UKBB/ldby_segment/snp_group2_maf1_chr",chr,".txt"), row.names=F, quote=F, col.names=F)
	write.table(lb2_snp_maf2, paste0("/scratch/users/tcheandj/subset.ukbb.for.clumping/ldscore_UKBB/ldby_segment/snp_group2_maf2_chr",chr,".txt"), row.names=F, quote=F, col.names=F)
	write.table(lb2_snp_maf3, paste0("/scratch/users/tcheandj/subset.ukbb.for.clumping/ldscore_UKBB/ldby_segment/snp_group2_maf3_chr",chr,".txt"), row.names=F, quote=F, col.names=F)
	write.table(lb2_snp_maf4, paste0("/scratch/users/tcheandj/subset.ukbb.for.clumping/ldscore_UKBB/ldby_segment/snp_group2_maf4_chr",chr,".txt"), row.names=F, quote=F, col.names=F)
	write.table(lb2_snp_maf5, paste0("/scratch/users/tcheandj/subset.ukbb.for.clumping/ldscore_UKBB/ldby_segment/snp_group2_maf5_chr",chr,".txt"), row.names=F, quote=F, col.names=F)


	#write.table(lb3_snp, paste0("/scratch/users/tcheandj/subset.ukbb.for.clumping/ldscore_UKBB/ldby_segment/snp_group3_chr",chr,".txt"), row.names=F, quote=F, col.names=F)
	write.table(lb3_snp_maf1, paste0("/scratch/users/tcheandj/subset.ukbb.for.clumping/ldscore_UKBB/ldby_segment/snp_group3_maf1_chr",chr,".txt"), row.names=F, quote=F, col.names=F)
	write.table(lb3_snp_maf2, paste0("/scratch/users/tcheandj/subset.ukbb.for.clumping/ldscore_UKBB/ldby_segment/snp_group3_maf2_chr",chr,".txt"), row.names=F, quote=F, col.names=F)
	write.table(lb3_snp_maf3, paste0("/scratch/users/tcheandj/subset.ukbb.for.clumping/ldscore_UKBB/ldby_segment/snp_group3_maf3_chr",chr,".txt"), row.names=F, quote=F, col.names=F)
	write.table(lb3_snp_maf4, paste0("/scratch/users/tcheandj/subset.ukbb.for.clumping/ldscore_UKBB/ldby_segment/snp_group3_maf4_chr",chr,".txt"), row.names=F, quote=F, col.names=F)
	write.table(lb3_snp_maf5, paste0("/scratch/users/tcheandj/subset.ukbb.for.clumping/ldscore_UKBB/ldby_segment/snp_group3_maf5_chr",chr,".txt"), row.names=F, quote=F, col.names=F)


	#write.table(lb4_snp, paste0("/scratch/users/tcheandj/subset.ukbb.for.clumping/ldscore_UKBB/ldby_segment/snp_group4_chr",chr,".txt"), row.names=F, quote=F, col.names=F)
	write.table(lb4_snp_maf1, paste0("/scratch/users/tcheandj/subset.ukbb.for.clumping/ldscore_UKBB/ldby_segment/snp_group4_maf1_chr",chr,".txt"), row.names=F, quote=F, col.names=F)
	write.table(lb4_snp_maf2, paste0("/scratch/users/tcheandj/subset.ukbb.for.clumping/ldscore_UKBB/ldby_segment/snp_group4_maf2_chr",chr,".txt"), row.names=F, quote=F, col.names=F)
	write.table(lb4_snp_maf3, paste0("/scratch/users/tcheandj/subset.ukbb.for.clumping/ldscore_UKBB/ldby_segment/snp_group4_maf3_chr",chr,".txt"), row.names=F, quote=F, col.names=F)
	write.table(lb4_snp_maf4, paste0("/scratch/users/tcheandj/subset.ukbb.for.clumping/ldscore_UKBB/ldby_segment/snp_group4_maf4_chr",chr,".txt"), row.names=F, quote=F, col.names=F)
	write.table(lb4_snp_maf5, paste0("/scratch/users/tcheandj/subset.ukbb.for.clumping/ldscore_UKBB/ldby_segment/snp_group4_maf5_chr",chr,".txt"), row.names=F, quote=F, col.names=F)

}

 sbatch -p jpriest --qos=normal -J snp_split -t 36:00:00 -N 1 -n 1 -o snp_split.log -e snp_split.err --mem=42000 R CMD BATCH  split.SNP.R 




