#!/bin/bash
#SBATCH -J list8  #Slurm job name
##SBATCH --mail-user=songdong_bio@126.com 
##SBATCH --mail-type=end
#SBATCH -p cpu-share
#SBATCH -o tra_P2.out
#SBATCH -e tra_P2.err
##SBATCH -N 1 -n 40
#SBATCH --cpus-per-task=40
#SBATCH --exclusive


#set -e
#source ~/.bashrc
#source ../savi.env

sh trans.mpi.sh P5 results_T_DNA >tmp8 &
wait
