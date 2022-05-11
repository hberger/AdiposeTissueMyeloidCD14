#!/bin/bash
#
##SBATCH --job-name=velocyto_BCN
#SBATCH --output=slurm_output.txt
#
##SBATCH --ntasks=12
#SBATCH --nodes=1
#SBATCH --time=60:00:00
#SBATCH --mem=80G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hilmar.berger@charite.de


SAMPLES=$(ls ../../../../cellranger/BCN_MarColl/Results)

for s in $SAMPLES; do

  INPUT_FOLDER=/data/BCN_MarColl/Results/$s
  GTF=/data/download/refdata-gex-GRCh38-2020-A/genes/genes.gtf

  singularity run --bind ~/work/cellranger:/data ~/work/velocyto.sif run10x -vv $INPUT_FOLDER $GTF

done
