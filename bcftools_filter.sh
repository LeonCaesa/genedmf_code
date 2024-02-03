#!/bin/bash -l 
# Cores 
#$ -pe omp 16  
#$ -binding linear:16   
#$ -j y  
#$ -P dmfgrp 
echo "=========================================================="  
echo "Starting on       : $(date)" 
echo "Running on node   : $(hostname)" 
echo "Current directory : $(pwd)"  
echo "Current job ID    : $JOB_ID" 
echo "Current job name  : $JOB_NAME"  
echo "Task index number : $TASK_ID" 
echo "==========================================================" 

module load bcftools
module load plink/2.0

#chr=19
chr=$SGE_TASK_ID
file="/projectnb2/dmfgrp/GeneDMF/data/chr${chr}.vcf.gz"
out1="/projectnb2/dmfgrp/GeneDMF/data/chr${chr}/chr${chr}_bcftools.vcf.gz"
out2="/projectnb2/dmfgrp/GeneDMF/data/chr${chr}/chr${chr}_pruned.vcf.gz"
mkdir /projectnb2/dmfgrp/GeneDMF/data/chr${chr}
cd /projectnb2/dmfgrp/GeneDMF/data/chr${chr}

bcftools view --types snps -q 0.01:minor -i 'F_MISSING<0.1' $file -Oz -o $out1

bcftools index $out1

plink2 --vcf $out1 --indep-pairwise 500 'kb' 1 0.1 --set-missing-var-ids @:#:\$r:\$a  --new-id-max-allele-len 2 missing --out chr${chr}

awk -F ':' '{print $1, $2}' OFS='\t'  chr${chr}.prune.in >> chr${chr}.prune_vars

bcftools view -R chr${chr}.prune_vars $out1 -Oz -o $out2

echo "=========================================================="   
echo "Finished on       : $(date)"   
echo "==========================================================" 
 
