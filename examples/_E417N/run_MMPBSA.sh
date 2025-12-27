#!/bin/bash
#SBATCH --job-name=run_66_mut
#SBATCH --partition=cpu
#SBATCH --output=run_mmpbsa_66.out
#SBATCH --error=run_mmpbsa_66.error
#SBATCH --time=48:00:00
echo "Loading modules..."
module load amber 
source amber22/amber.sh
$AMBERHOME/bin/MMPBSA.py -O -i mmpbsa.in -o FINAL_RESULTS_MMPBSA_6m0j_noHetE417N.dat -sp 6m0j_noHetE417N_solvated.prmtop -cp 6m0j_noHetE417N.prmtop -rp 6m0j_noHetE417N_recpt.prmtop -lp 6m0j_noHetE417N_ligand.prmtop -y *.mdcrd
