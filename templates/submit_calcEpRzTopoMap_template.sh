#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=username@mail.com
#SBATCH --export=ALL
#SBATCH -J calcEpRzTopoMap_template
#SBATCH --partition=xeon24
#SBATCH -n 24
#SBATCH -N 1
#SBATCH --mem=200G
#SBATCH --time=08:00:00
#SBATCH --output=log_file_calcEpRzTopoMap_template.out

julia start_calcEpRzTopoMap_template.jl
mv start_calcEpRzTopoMap_template.jl ../OWCF_results/template/
mv log_file_calcEpRzTopoMap_template.out ../OWCF_results/template/
mv submit_calcEpRzTopoMap_template.sh ../OWCF_results/template/
