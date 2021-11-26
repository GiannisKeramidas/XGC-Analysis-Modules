#!/bin/bash
#SBATCH -A mp118
#SBATCH -N 60
#SBATCH -n 60
#SBATCH -p regular
#SBATCH -J par_job
#SBATCH --mail-user=giannis.kx@gmail.com
#SBATCH -t 2:00:00
module load python
module load mpi4py
srun -n 60 python /global/homes/g/giannos/xgc_python_dir/parallel.py exb_flux
