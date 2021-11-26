# setup script for accessing XGC1 output via adios
# usage:
#	source setupad.sh
#

# to refresh xgc.py:
#	cd python_xgc
#	git pull


#module load python # loads python_2.7.9 and numpy,scipy,matplotlib
#module load python/2.7-anaconda # Michael recommends (at least for h5)
#module unload python # if you want to run setup.sh

# for adios
#module use -a /project/projectdirs/m499/jyc/edison/sw/modulefiles
module use -a /global/homes/g/giannos/modulefiles
module load python
module load adios/1.8.1-pre
module load python_adios

PYTHONPATH='/global/homes/g/giannos/xgc_python_dir:'$PYTHONPATH
export PYTHONPATH
