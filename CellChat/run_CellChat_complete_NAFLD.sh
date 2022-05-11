#!/bin/bash
#
#SBATCH --job-name=CellChat_complete_NAFLD
#SBATCH --output=slurm_output_complete_NAFLD
#
#SBATCH --ntasks=12
#SBATCH --nodes=1
#SBATCH --time=80:00:00
#SBATCH --mem=80G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hilmar.berger@charite.de

r_file=run_CellChat_complete_NAFLD.R

(cat <<"EOF"
library(knitr)
rmarkdown::render("CellChat_complete_NAFLD.Rmd")
quit("no")
EOF
) > $r_file


R CMD BATCH $r_file

