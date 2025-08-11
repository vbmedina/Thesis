#PBS -lwalltime=36:00:00
#PBS -lselect=1:ncpus=16:mem=128gb:ngpus=1
#PBS -o /rds/general/user/vbm24/home/Thesis/p2_models/models/outputs
#PBS -e /rds/general/user/vbm24/home/Thesis/p2_models/models/outputs
#PBS -N dmpnn_norm_agg
 
 
source ${HOME}/.bashrc
conda activate chemprop
 
cd /rds/general/user/qg622/home/Y3/chemprop/chemprop
 
cell_names=(
    '786-0' 'A498' 'A549_ATCC' 'ACHN' 'BT-549' 'CAKI-1'
    'CCRF-CEM' 'COLO_205' 'DU-145' 'EKVX' 'HCC-2998' 'HCT-116'
    'HCT-15' 'HL-60_TB_' 'HOP-62' 'HOP-92' 'HS_578T' 'HT29'
    'IGROV1' 'K-562' 'KM12' 'LOX_IMVI' 'M14' 'MALME-3M'
    'MCF7' 'MDA-MB-231_ATCC' 'MDA-MB-435' 'MDA-N' 'MOLT-4'
    'NCI-H226' 'NCI-H23' 'NCI-H322M' 'NCI-H460' 'NCI-H522'
    'NCI_ADR-RES' 'OVCAR-3' 'OVCAR-4' 'OVCAR-5' 'OVCAR-8'
    'PC-3' 'RPMI-8226' 'RXF_393' 'SF-268' 'SF-295' 'SF-539'
    'SK-MEL-2' 'SK-MEL-28' 'SK-MEL-5' 'SK-OV-3' 'SN12C'
    'SNB-19' 'SNB-75' 'SR' 'SW-620' 'T-47D' 'TK-10' 'U251'
    'UACC-257' 'UACC-62' 'UO-31'
)
 
 
for seed in 42 43 44 45 46; do
    for cell_line in "${cell_names[@]}"; do
        for split in umap random scaffold butina; do
            predictions_file="/rds/general/user/qg622/home/Y3/nature/results/D-MPNN_norm_agg/$cell_line/split_${split}/regression_$seed/test_predictions.csv"
            if [[ ! -f "$predictions_file" ]]; then
                python run_model.py --cell_line "$cell_line" --seed "$seed" --split "$split"
                # python run_model.py --cell_line "$cell_line" --seed "$seed" --split "umap"
            else
                echo "Predictions file already exists for $cell_line with seed $seed. Skipping..."
            fi
        done
    done
done