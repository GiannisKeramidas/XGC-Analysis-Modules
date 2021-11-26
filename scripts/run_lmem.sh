#!/bin/bash
#SBATCH -A m499
#SBATCH --clusters=escori
#SBATCH --qos=bigmem
#SBATCH --nodes=1
#SBATCH --time=08:00:00
#SBATCH --job-name=big_memory
#SBATCH --mem=500GB

srun -n 1 ./large_mem.py > lmem.out
