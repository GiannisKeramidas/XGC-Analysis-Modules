#!/bin/bash
#SBATCH -A mp118
#SBATCH -N 40
#SBATCH -C haswell
#SBATCH -n 80
#SBATCH -p debug
#SBATCH -J par_job
#SBATCH --mail-user=giannis.kx@gmail.com
#SBATCH -t 00:30:00
module load python/2.7-anaconda
srun -n 80 -c 32 --cpu_bind=cores python /global/homes/g/giannos/xgc_python_dir/parallel.py shear_calc
