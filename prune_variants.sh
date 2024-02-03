#!/bin/bash -l 
#$ -l h_vmem=10G 
#$ -l h_rt=12:00:00
#$ -l cpu_type=X5670
# -l mem_total=93G
# -l mem_per_core=8G
# Cores 
#$ -pe smp 4
#$ -binding linear:4   
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

#cd /rprojectnb2/necs/Sophia_Analysis/PRS_repo/LLFS 
OMP_NUM_THREADS=1

module load plink/2.0
module load R 

#chr=$SGE_TASK_ID
chr=22
file="/projectnb2/dmfgrp/GeneDMF/data/chr${chr}.vcf.gz"
out="/projectnb2/dmfgrp/GeneDMF/data/chr${chr}/chr${chr}_filtered"

if [ -d "/projectnb2/dmfgrp/GeneDMF/data/chr${chr}/"]; then
  rm -rf /projectnb2/dmfgrp/GeneDMF/data/chr${chr}/
fi
mkdir /projectnb2/dmfgrp/GeneDMF/data/chr${chr}

plink --vcf $file --indep-pairwise 500 'kb' 1 0.1 --set-missing-var-ids @_#  --new-id-max-allele-len 2 missing  --out $out 
#--threads 4
echo "=========================================================="   
echo "Finished on       : $(date)"   
echo "==========================================================" 







