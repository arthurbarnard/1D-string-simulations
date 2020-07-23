#!/bin/bash
#SBATCH --job-name=CNT_test
#SBATCH --account=genpool
#SBATCH --partition=hpc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=3:00:00
#SBATCH --mem=1G
#SBATCH --chdir=/gscratch/scrubbed/wbhuynh/1D-string-simulations
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wbhuynh@uw.edu
#SBATCH --export=all
./CNT_run.exe 10 ./initfiles/inputFile testOut ./initfiles/seedTable.txt
