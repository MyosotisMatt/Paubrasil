#!/bin/bash
#SBATCH --job-name=Admixture5
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=edeline.gagnon@gmail.com
#SBATCH --partition=long


#Removes multiallelic alleles
#bcftools view --max-alleles 2 -o data5moda.vcf --exclude-types indels data5mod.vcf

cd /mnt/shared/scratch/egagnon/PauBrasil/07_Admixture/data5mod3/

## Pruning VCF file for admixture analysis. (This might not be required for you @Edeline)
plink --vcf data5moda.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 10 0.5 --out prunned_data5moda

#This is similar, but makes something that can then be used for PCA
plink --vcf data5moda.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# \
--extract prunned_data5moda.prune.in \
--make-bed --pca --out paubrasil

eval "$(conda shell.bash hook)"
conda activate admixture

# make bed file with pruned data (Start from here with your VCF file)
plink --vcf data5moda.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --extract prunned_data5moda.prune.in --make-bed --out prunned_data

# ADMIXTURE does not accept chromosome names that are not human chromosomes. Exchange the first column by 0

awk '{$1="0";print $0}' prunned_data.bim > prunned_data.bim.tmp
mv prunned_data.bim.tmp prunned_data.bim

# You can run this code which is the standard one.
for K in {1..10} 
do 
admixture --cv=10 prunned_data.bed $K -j2 | tee log${K}.out 
done

# If you want further repetition of each K, then use this code.
for K in {1..10} 
do 
	for repeat in {1..10} 
	do 
		admixture --cv prunned_data.bed $K -j16 -s ${repeat} | tee log${K}.${repeat}.out
	done
done

FILE=prunned_data

awk '/CV/ {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' *out | cut -c 4,7-20 > prunned_data.cv.error
awk '{split($1,name,"."); print $1,name[2]}' ${FILE}.nosex > prunned_data.list

#Rscript admixture_script.R -p Rscript_results -i prunned_data.list -k 10 -l Paubrasilia_echinata_22_7905_13,Paubrasilia_echinata_28_7905_10,Paubrasilia_echinata_29_7905_06,Paubrasilia_echinata_31_7905_04,Paubrasilia_echinata_32_7905_07,Paubrasilia_echinata_33_7905_03,Paubrasilia_echinata_47_7905_14,Paubrasilia_echinata_43_7905_01,Paubrasilia_echinata_34_7905_02,Paubrasilia_echinata_30_7905_05,Paubrasilia_echinata_65_7905_11,Paubrasilia_echinata_23_7905_12,Paubrasilia_echinata_49_472149,Paubrasilia_echinata_69_7910_09,Paubrasilia_echinata_52_7596,Paubrasilia_echinata_57_7595,Paubrasilia_echinata_64_7987_02,Paubrasilia_echinata_60_7893_01,Paubrasilia_echinata_61_7894_04,Paubrasilia_echinata_62_7896_03,Paubrasilia_echinata_63_7898_05,Paubrasilia_echinata_41_1624,Paubrasilia_echinata_27_7923_01,Paubrasilia_echinata_25_7924_03,Paubrasilia_echinata_26_7923_02,Paubrasilia_echinata_45_7924_04,Paubrasilia_echinata_42_6659,Paubrasilia_echinata_45_7918_02,Paubrasilia_echinata_66_7918_03,Paubrasilia_echinata_70_7918_01,Paubrasilia_echinata_59_6873,Paubrasilia_echinata_67_7918_04,Paubrasilia_echinata_68_7918_05,Paubrasilia_echinata_56_6662,Paubrasilia_echinata_51_6961,Paubrasilia_echinata_58_365034,Paubrasilia_echinata_53_12IX1985,Paubrasilia_echinata_38_Pereira_365015,Paubrasilia_echinata_40_8123
#Rscript admixture_script.R -p Rscript_results -i prunned_data.list -k 10 -l Paubrasilia_echinata_22,Paubrasilia_echinata_23_7905_12,Paubrasilia_echinata_25_7924_03,Paubrasilia_echinata_26_7923_02,Paubrasilia_echinata_27_7923_01,Paubrasilia_echinata_28_7905_10,Paubrasilia_echinata_29_7905_06,Paubrasilia_echinata_30_7905_05,Paubrasilia_echinata_31_7905_04,Paubrasilia_echinata_32_7905_07,Paubrasilia_echinata_33_7905_03,Paubrasilia_echinata_34_7905_02,Paubrasilia_echinata_38_Pereira_365015,Paubrasilia_echinata_40_8123,Paubrasilia_echinata_41_1624,Paubrasilia_echinata_42_6659,Paubrasilia_echinata_43_7905_01,Paubrasilia_echinata_45_7924_04,Paubrasilia_echinata_46_7918_02,Paubrasilia_echinata_47_7905_14,Paubrasilia_echinata_49_472149,Paubrasilia_echinata_51_6961,Paubrasilia_echinata_52_7596,Paubrasilia_echinata_53_12IX1985,Paubrasilia_echinata_56_6662,Paubrasilia_echinata_57_7595,Paubrasilia_echinata_58_365034,Paubrasilia_echinata_59_6873,Paubrasilia_echinata_60_7893_01,Paubrasilia_echinata_61_7894_04,Paubrasilia_echinata_62_7896_03,Paubrasilia_echinata_63_7898_05,Paubrasilia_echinata_64_7987_02,Paubrasilia_echinata_65_7905_11,Paubrasilia_echinata_66_7918_03,Paubrasilia_echinata_67_7918_04,Paubrasilia_echinata_68_7918_05,Paubrasilia_echinata_69_7910_09,Paubrasilia_echinata_70_7918_01

Rscript admixture_script.R -p prunned_data -i prunned_data.list -o repeat_results_prunned_data5mod3 -k 10 -l 7905,7910,472149,7893,7894,7896,7898,7987,7595,7596,7924,7923,1624,Pereira365015,8123,6961,12IX1985,6662,365034,6873,6659,7918

#awk '/CV/ {print $5}' *out | cut -c 4,7-20 > prunned_data.cv.error
#awk '{split($1,name,"."); print $1,name[2]}' ${FILE}.nosex > prunned_data.list

#Rscript admixture_script.R -p prunned_data -i prunned_data.list -o repeat_results_prunned_data5moda_k5 -k 10 -l 7905,7910,472149,7893,7894,7896,7898,7987,7595,7596,7924,7923,1624,Pereira365015,8123,6961,12IX1985,6662,365034,6873,6659,7918


########################
#
#plink --vcf data5moda.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --extract prunned_data5moda.prune.in --out data5modeb

#Transform recode bed file into vcf
#plink --bfile prunned_data --recode vcf --out prunned_data_mode_v2

