#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=username@mail.com
#SBATCH --export=ALL
#SBATCH -J topoMap_template
#SBATCH --partition=xeon24
#SBATCH -n 24
#SBATCH -N 1
#SBATCH --mem=150G
#SBATCH --time=04:00:00
#SBATCH --output=log_file_calcTopoMap_template.out
julia start_calcTopoMap_template.jl
mv start_calcTopoMap_template.jl ../OWCF_results/template
mv log_file_calcTopoMap_template.out ../OWCF_results/template
mv submit_calcTopoMap_template.sh ../OWCF_results/template
