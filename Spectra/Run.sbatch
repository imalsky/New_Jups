#!/bin/bash
#SBATCH --account=pi-larogers
#SBATCH --job-name=MalskyTest
#SBATCH --time=36:00:00
##SBATCH --partition=sandyb
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=isaacmalsky@gmail.com

module load python/anaconda-2020.02
source activate /home/aneil/conda
module unload python

export OMP_NUM_THREADS=1

python run_spectra.py