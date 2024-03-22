#! /bin/bash
#$ -S /bin/bash
#$ -q NOParalela
#$ -l NOParalela
#$ -t 1:15
module load alhambra/R-3.5.2
R CMD BATCH Type_II_error.R Type_II_error.Rout.$SGE_TASK_ID