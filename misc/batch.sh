#!/bin/bash

#SBATCH -o /home/hpc/pr63so/ga25cux2/turbulent/script_output.%j.out
#SBATCH -D /home/hpc/pr63so/ga25cux2/
#SBATCH -J turbm1.5
#SBATCH --mail-type=END
#SBATCH --mail-user=johannes.klicpera@tum.de
#SBATCH --export=NONE
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --partition=snb
#SBATCH --begin=now

source /etc/profile.d/modules.sh
module load python

export OMP_NUM_THREADS=16

cd /home/hpc/pr63so/ga25cux2/turbulent

mpirun -np 1 ./ns conf/channel/m1.5.xml
