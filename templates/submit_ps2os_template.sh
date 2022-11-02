#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=username@mail.com
#SBATCH --export=ALL
#SBATCH -J ps2os_template
#SBATCH --partition=xeon24
#SBATCH -n 24
#SBATCH -N 1
#SBATCH --mem=250G
#SBATCH --time=06:00:00
#SBATCH --output=log_file_ps2os_template.out

julia start_ps2os_template.jl
mv start_ps2os_template.jl ../OWCF_results/template/
mv log_file_ps2os_template.out ../OWCF_results/template/
mv submit_ps2os_template.sh ../OWCF_results/template/
