#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=username@mail.com
#SBATCH --export=ALL
#SBATCH -J calcOG_template
#SBATCH --partition=xeon24
#SBATCH -n 24
#SBATCH -N 1
#SBATCH --mem=100G
#SBATCH --time=01:00:00
#SBATCH --output=log_file_calcOG_template.out

julia start_calcOG_template.jl
mv start_calcOG_template.jl ../OWCF_results/template/
mv log_file_calcOG_template.out ../OWCF_results/template/
mv submit_calcOG_template.sh ../OWCF_results/template/
