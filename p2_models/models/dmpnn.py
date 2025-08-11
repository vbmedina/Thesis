"""
scripts/train_dmpnn.py

Train a Chemprop D-MPNN on your pp_tidy.csv with 5-fold cross-validation.
"""
import sys
import chemprop
import args
from chemprop.args import TrainArgs
from chemprop.train import cross_validate, run_training

def main():
    # Base arguments — tweak as needed
    base_args = [
        '--data_path',    'data/tidy/pp_tidy.csv',
        '--dataset_type', 'regression',           # or classification, multiclass, spectra
        '--save_dir',     'checkpoints/dmpnn',
        '--split_type',   'scaffold_balanced',    # random, random_stratified, scaffold_balanced, fingerprint, butina, umap, etc.
        '--split_sizes',  '0.8', '0.1', '0.1',
        '--num_folds',    '5',
        '--seed',         '42',
        '--features_generator', 'rdkit_2d_normalized',
        '--cache_cutoff', 'inf',                  # keep all graphs in memory
        '--quiet'                                 # remove to see full logs
    ]

    # Allow CLI overrides
    args = TrainArgs().parse_args(base_args + sys.argv[1:])

    mean_score, std_score = cross_validate(
        args=args,
        train_func=run_training
    )
    print(f'5-fold CV result: {mean_score:.4f} ± {std_score:.4f}')

if __name__ == '__main__':
    main()