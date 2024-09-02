#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=username@mail.com
#SBATCH --export=ALL
#SBATCH -J calc4DW_template
#SBATCH --partition=xeon24
#SBATCH -n 24
#SBATCH -N 1
#SBATCH --mem=250G
#SBATCH --time=2-00:00:00
#SBATCH --output=log_file_calc4DW_template.out

julia start_calc4DW_template.jl
mv start_calc4DW_template.jl ../OWCF_results/template/
mv log_file_calc4DW_template.out ../OWCF_results/template/
mv submit_calc4DW_template.sh ../OWCF_results/template/
