#!/bin/bash
#SBATCH --job-name=CNT_test
#SBATCH --output=CNT_test.%j.out
#SBATCH --error=CNT_test.%j.out
#SBATCH --account=genpool
#SBATCH --partition=hpc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=3
#SBATCH --time=48:00:00
#SBATCH --mem=1G
#SBATCH --chdir=/gscratch/scrubbed/wbhuynh/1D-string-simulations
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=wbhuynh@uw.edu
#SBATCH --export=all

((jobNumIn=$SLURM_ARRAY_TASK_ID))
((listNum = jobNumIn*3))

parallel --delay .2 -j $SLURM_NTASKS ./CNT_run.exe {1} ./initfiles/inputFile testOut ./initfiles/seedTable.txt ./initfiles/mass ./initfiles1/K_cnt ./initfiles/kap_cnt ::: {$listNum,$(($listNum+1)),$(($listNum+2))}
