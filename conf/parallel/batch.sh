#!/bin/bash

#SBATCH -o /home/hpc/pr63so/ga25cux2/turbulent/script_output.%j.out
#SBATCH -D /home/hpc/pr63so/ga25cux2/
#SBATCH -J turbulent
#SBATCH --export=NONE
#SBATCH --time=2:00:00
#SBATCH --nodes={nodes}
#SBATCH --partition=snb
#SBATCH --constraint=turbo_off
#SBATCH --begin=now

source /etc/profile.d/modules.sh

export OMP_NUM_THREADS=16

cd /home/hpc/pr63so/ga25cux2/turbulent
