#!/bin/bash
#
#SBATCH --job-name=CPDB_BCN
#SBATCH --output=slurm_output.txt
#
#SBATCH --ntasks=12
#SBATCH --nodes=1
#SBATCH --time=80:00:00
#SBATCH --mem=80G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hilmar.berger@charite.de

DATA_FOLDER=~/work/data/scCell/BCN_2021/Analysis/Results/data_storage/2021-05-31

cellphonedb method statistical_analysis $DATA_FOLDER/BCN_MoMF_DC_only_meta_CPDB.txt $DATA_FOLDER/BCN_MoMF_DC_only_norm_counts_CPDB.txt --counts-data hgnc_symbol

cellphonedb plot dot_plot

cellphonedb plot heatmap_plot $DATA_FOLDER/BCN_MoMF_DC_only_meta_CPDB.txt
