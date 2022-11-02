#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=username@mail.com
#SBATCH --export=ALL
#SBATCH -J extrNullOrbs_template
#SBATCH --partition=xeon24
#SBATCH -n 24
#SBATCH -N 1
#SBATCH --mem=150G
#SBATCH --time=02:00:00
#SBATCH --output=log_file_extrNullOrbs_template.out
julia start_extrNullOrbs_template.jl
mv start_extrNullOrbs_template.jl ../OWCF_results/template
mv log_file_extrNullOrbs_template.out ../OWCF_results/template
mv submit_extrNullOrbs_template.sh ../OWCF_results/template
