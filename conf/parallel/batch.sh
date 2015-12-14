#!/bin/bash

#SBATCH -o /home/hpc/pr63so/ga25cux2/turbulent/output/{run_name}/script_output.%j.out
#SBATCH -D /home/hpc/pr63so/ga25cux2/
#SBATCH -J {run_name}
#SBATCH --export=NONE
#SBATCH --time=1:00:00
#SBATCH --nodes={nodes}
#SBATCH --partition=mpp2_inter
# #SBATCH --constraint=turbo_off
#SBATCH --begin=now

source /etc/profile.d/modules.sh

export OMP_NUM_THREADS=16

cd /home/hpc/pr63so/ga25cux2/turbulent
