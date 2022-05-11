#!/bin/bash
#
#SBATCH --job-name=CellChat_diseased_Ramachandran_and_NAFLD_adipose_myeloid_cells_secreted_ligands_v2
#SBATCH --output=slurm_output_CellChat_diseased_Ramachandran_and_NAFLD_adipose_myeloid_cells_secreted_ligands_v2
#
#SBATCH --ntasks=12
#SBATCH --nodes=1
#SBATCH --time=80:00:00
#SBATCH --mem=80G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hilmar.berger@charite.de

r_file=run_CellChat_diseased_Ramachandran_and_NAFLD_adipose_myeloid_cells_secreted_ligands_v2.R

(cat <<"EOF"
library(knitr)
rmarkdown::render("CellChat_NAFLD_with_Ramachandran_secreted_ligands_cirrhotic_liver_and_NAFLD_adipose_tissue_only_v2.Rmd")
quit("no")
EOF
) > $r_file


R CMD BATCH $r_file

