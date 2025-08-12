#PBS -lwalltime=36:00:00
#PBS -lselect=1:ncpus=16:mem=128gb:ngpus=1
#PBS -o /rds/general/user/vbm24/home/Thesis/p2_models/models/outputs
#PBS -e /rds/general/user/vbm24/home/Thesis/p2_models/models/outputs
#PBS -N dmpnn_norm_agg
 
 
source ${HOME}/.bashrc
conda activate chemprop
 
cd /rds/general/user/vbm24/home/Thesis/p2_models/models

python training.py