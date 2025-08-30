#PBS -lwalltime=36:00:00
#PBS -lselect=1:ncpus=16:mem=128gb:ngpus=1
#PBS -o /rds/general/user/vbm24/home/Thesis/p2_models/models/outputs
#PBS -e /rds/general/user/vbm24/home/Thesis/p2_models/models/outputs
#PBS -N dmpnn_train_val
 
 
source ${HOME}/.bashrc
conda activate /rds/general/user/vbm24/home/Thesis/.conda/dmpnn
 
cd "/rds/general/user/vbm24/home/Thesis/p2_models/1 - dmpnn_train_val_test_code"

python "./1 - train_validate_dmpnn.py"