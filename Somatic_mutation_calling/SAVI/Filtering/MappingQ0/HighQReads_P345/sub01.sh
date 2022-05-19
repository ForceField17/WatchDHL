#!/bin/bash
#SBATCH -J split01  #Slurm job name
##SBATCH --mail-user=songdong_bio@126.com 
##SBATCH --mail-type=end
#SBATCH -p cpu-share
#SBATCH -o tmp01.out
#SBATCH -e tmp01.err
##SBATCH -N 1 -n 40
#SBATCH --cpus-per-task=40
#SBATCH --exclusive


#set -e
#source ~/.bashrc
source ../../../../../savi.env
sh trans.mpi.sh split01 results_T_DNA
wait
