#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=username@mail.com
#SBATCH --export=ALL
#SBATCH -J 4Dto2D_template
#SBATCH --partition=xeon24
#SBATCH -n 24
#SBATCH -N 1
#SBATCH --mem=100G
#SBATCH --time=02:00:00
#SBATCH --output=log_file_4Dto2D_template.out

julia start_4Dto2D_template.jl
mv start_4Dto2D_template.jl ../OWCF_results/template/
mv log_file_4Dto2D_template.out ../OWCF_results/template/
mv submit_4Dto2D_template.sh ../OWCF_results/template/
