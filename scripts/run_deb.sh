#!/bin/bash
#SBATCH -A mp118
#SBATCH -N 60
#SBATCH -C haswell
#SBATCH -p debug
#SBATCH -J par_job
#SBATCH --mail-user=giannis.kx@gmail.com
#SBATCH -t 00:30:00
module load python/2.7-anaconda
srun -n 60 -c64 --cpu_bind=cores python /global/homes/g/giannos/xgc_python_dir/parallel.py
