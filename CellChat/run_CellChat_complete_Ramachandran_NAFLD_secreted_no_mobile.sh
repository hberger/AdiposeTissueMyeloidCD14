#!/bin/bash
#
#SBATCH --job-name=CellChat_complete_Ramachandran_NAFLD_secreted_no_mobile
#SBATCH --output=slurm_output_complete_Ramachandran_NAFLD_secreted_no_mobile
#
#SBATCH --ntasks=12
#SBATCH --nodes=1
#SBATCH --time=80:00:00
#SBATCH --mem=80G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hilmar.berger@charite.de

r_file=run_CellChat_complete_Ramachandran_NAFLD_secreted_no_mobile.R

(cat <<"EOF"
library(knitr)
rmarkdown::render("CellChat_NAFLD_with_Ramachandran_secreted_ligands_residential_populations_only.Rmd")
quit("no")
EOF
) > $r_file


R CMD BATCH $r_file

