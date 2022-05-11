#!/bin/bash
#
#SBATCH --job-name=CellChat_healthy_Ramachandran_and_healthy_adipose_myeloid_cells_secreted_ligands
#SBATCH --output=slurm_output_CellChat_NAFLD_with_Ramachandran_secreted_ligands_healthy_liver_and_adipose_tissue_only.Rmd
#
#SBATCH --ntasks=12
#SBATCH --nodes=1
#SBATCH --time=80:00:00
#SBATCH --mem=80G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hilmar.berger@charite.de

r_file=run_CellChat_healthy_Ramachandran_and_healthy_adipose_myeloid_cells_secreted_ligands.R

(cat <<"EOF"
library(knitr)
rmarkdown::render("CellChat_NAFLD_with_Ramachandran_secreted_ligands_healthy_liver_and_adipose_tissue_only.Rmd")
quit("no")
EOF
) > $r_file


R CMD BATCH $r_file

