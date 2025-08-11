#!/bin/bash
#PBS -l select=1:ncpus=2:mem=4gb
#PBS -l walltime=00:01:00
#PBS -N hello_world

module load Python/3.12.3-GCCcore-13.3.0

cd $PBS_O_WORKDIR

python /rds/general/user/vbm24/home/Thesis/p2_models/models/dmpnn.py