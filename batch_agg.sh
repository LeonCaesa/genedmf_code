#!/bin/bash -l 
# Cores 
#$ -j y  
#$ -P dmfgrp
#$ -l h_rt=72:00:00

echo "=========================================================="  
echo "Starting on       : $(date)" 
echo "Running on node   : $(hostname)" 
echo "Current directory : $(pwd)"  
echo "Current job ID    : $JOB_ID" 
echo "Current job name  : $JOB_NAME"  
echo "Task index number : $TASK_ID" 
echo "=========================================================="   


module load R/4.0.2
Rscript /projectnb/dmfgrp/GeneDMF/agg_script.R


echo "=========================================================="   
echo "Finished on       : $(date)"   
echo "==========================================================" 







