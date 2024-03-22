#! /bin/bash
#$ -S /bin/bash
#$ -q NOParalela
#$ -l NOParalela
#$ -t 1:15
module load alhambra/R-3.5.2
R CMD BATCH Type_I_error.R Type_I_error.Rout.$SGE_TASK_ID