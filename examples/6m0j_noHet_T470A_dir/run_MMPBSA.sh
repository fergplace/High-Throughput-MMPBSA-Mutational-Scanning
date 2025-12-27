#!/bin/bash
#SBATCH --job-name=run_6m0j_noHet_mut
#SBATCH --partition=cpu
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem=10000
#SBATCH --output=run_mmpbsa_6m0j_noHet.out
#SBATCH --error=run_mmpbsa_6m0j_noHet.error
#SBATCH --time=72:00:00
echo "Loading modules..."
module load amber 
source amber22/amber.sh

tleap -s -f tleap_mut.in > tleap_mut.out
$AMBERHOME/bin/MMPBSA.py -O -i C:\Users\13108\Documents\GitHub\HTMS_Amber\examples/mmpbsa.in -o FINAL_RESULTS_MMPBSA_tleap_6m0j_noHet_T470A.dat -sp 6m0j_noHet_solvated.prmtop -cp 6m0j_noHet.prmtop -rp 6m0j_noHet_recpt.prmtop -lp 6m0j_noHet_ligand.prmtop -y C:\Users\13108\Documents\GitHub\HTMS_Amber\examples/*.mdcrd -mc 6m0j_noHet_T470A.prmtop -ml 6m0j_noHet_T470A_ligand.prmtop
