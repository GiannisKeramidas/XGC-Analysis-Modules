#!/bin/bash
#SBATCH -A mp118
#SBATCH -N 3
#SBATCH -C haswell
#SBATCH -n 80
#SBATCH -p regular
#SBATCH -J par_job
#SBATCH --mail-user=giannis.kx@gmail.com
#SBATCH -t 04:00:00
module load python/2.7-anaconda
srun -n 80 -c 2 --cpu_bind=cores python /global/homes/g/giannos/xgc_python_dir/parallel.py shear_calc

