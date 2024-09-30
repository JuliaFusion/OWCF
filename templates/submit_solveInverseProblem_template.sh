# THIS SCRIPT IS CURRENTLY UNDER CONSTRUCTION!!!
# THIS SCRIPT IS CURRENTLY UNDER CONSTRUCTION!!!
# THIS SCRIPT IS CURRENTLY UNDER CONSTRUCTION!!!
#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=username@mail.com
#SBATCH --export=ALL
#SBATCH -J solveInverseProblem_template
#SBATCH --partition=xeon24
#SBATCH -n 24
#SBATCH -N 1
#SBATCH --mem=250G
#SBATCH --time=2-00:00:00
#SBATCH --output=log_file_solveInverseProblem_template.out

julia start_solveInverseProblem_template.jl
mv start_solveInverseProblem_template.jl ../OWCF_results/template/
mv log_file_solveInverseProblem_template.out ../OWCF_results/template/
mv submit_solveInverseProblem_template.sh ../OWCF_results/template/
