#!/bin/bash


## args: seurat object name, column in metadata to use, output directory, comparisons_file
## comparisons_file should be space delimited text file of format: cell_type1 celltype_2
## each pairwise comparison of interest should occur on separate line. No directionality- i.e., celltype1 --> celltype2  +  celltype2 --> celltype1 interactions will be tested
## downsample: TRUE/FALSE. Whether or not to downsample each program to the minimum amount of cells in a program from the same broad subtype
## if downsampling: broad subtype (tumor/tcell/myeloid) needs to be in metadata column of seurat objects as "highLevel_annot"

if [[ ${#} -ne 5 ]]
then 
        echo "usage: source CellChat-batch.sh seurat_obj CellAnnot_column_in_seurat output_dir comparisons_file downsample"
        #exit 1
fi 

seurat_obj=${1}
CellAnnot_Column=${2}
output_dir=${3}
comparison_file=${4}
downsample=${5}

script_dir=/home/jjl78/scripts/R/CellChat/

while read -r line; do 
	comparison=($line)
	celltype_1=${comparison[0]}
	celltype_2=${comparison[1]}
	name="${celltype_1}.vs.${celltype_2}"	
	echo "Running: ${celltype_1} vs ${celltype_2}"
	sbatch -o ${output_dir}/CC_${name}.out -e ${output_dir}/CC_${name}.err ${script_dir}/CellChat.sbatch ${seurat_obj} ${CellAnnot_Column} ${output_dir} ${celltype_1} ${celltype_2} ${downsample}

done < ${comparison_file}


