#!/bin/bash
#SBATCH -J split03  #Slurm job name
##SBATCH --mail-user=songdong_bio@126.com 
##SBATCH --mail-type=end
#SBATCH -p x-gpu-share
#SBATCH -o tmp03.out
#SBATCH -e tmp03.err
##SBATCH -N 1 -n 40
#SBATCH --cpus-per-task=20
#SBATCH --exclusive


#set -e
#source ~/.bashrc
source ../../../../../savi.env
sh trans.mpi.sh split03 results_T_DNA
wait
