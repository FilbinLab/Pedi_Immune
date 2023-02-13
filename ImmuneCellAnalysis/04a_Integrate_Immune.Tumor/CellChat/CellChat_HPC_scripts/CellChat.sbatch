#!/bin/bash

#SBATCH -p short 		# partition name
#SBATCH -t 0-10:00 		# hours:minutes runlimit after which job will be killed
#SBATCH --mem 12G 		# amount of memory requested
#SBATCH --job-name CellChat_R_script 		# Job name
#SBATCH -o %j.out		# File to which standard out will be written
#SBATCH -e %j.err 		# File to which standard err will be written

# Load required modules
module load gcc/6.2.0 R/4.0.1

# Point to personal library, if required
export R_LIBS_USER="~/R/4.0.1/library"


# Run the R script
Rscript ~/scripts/R/CellChat/CellChat.R ${1} ${2} ${3} ${4} ${5} ${6}

## args: seurat object name, column in metadata to use, output directory, cell type of interest #1, cell type of interest #2, whether or not to downsample (TRUE/FALSE)
