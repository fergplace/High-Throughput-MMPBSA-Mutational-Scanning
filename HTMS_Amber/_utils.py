import os
import sys

amino_acids = {
    'A': 'ALA',  # Alanine
    'R': 'ARG',  # Arginine
    'N': 'ASN',  # Asparagine
    'D': 'ASP',  # Aspartic Acid
    'C': 'CYS',  # Cysteine
    'E': 'GLU',  # Glutamic Acid
    'Q': 'GLN',  # Glutamine
    'G': 'GLY',  # Glycine
    'H': 'HIS',  # Histidine
    'I': 'ILE',  # Isoleucine
    'L': 'LEU',  # Leucine
    'K': 'LYS',  # Lysine
    'M': 'MET',  # Methionine
    'F': 'PHE',  # Phenylalanine
    'P': 'PRO',  # Proline
    'S': 'SER',  # Serine
    'T': 'THR',  # Threonine
    'W': 'TRP',  # Tryptophan
    'Y': 'TYR',  # Tyrosine
    'V': 'VAL',  # Valine
}
#inverse map 
three_letter_amino = {
    'ALA': 'A',
    'ARG': 'R',
    'ASN': 'N',
    'ASP': 'D',
    'CYS': 'C',
    'GLU': 'E',
    'GLN': 'Q',
    'GLY': 'G',
    'HIS': 'H',
    'ILE': 'I',
    'LEU': 'L',
    'LYS': 'K',
    'MET': 'M',
    'PHE': 'F',
    'PRO': 'P',
    'SER': 'S',
    'THR': 'T',
    'TRP': 'W',
    'TYR': 'Y',
    'VAL': 'V'
}


def get_amber_home_path():
    """
    Retrieves the AMBERHOME environment variable.
    Raises an error if AMBERHOME is not set.
    """
    amber_home = os.environ.get('AMBERHOME')
    if amber_home is None:
        raise EnvironmentError(
            "AMBERHOME environment variable is not set. "
            "Please ensure you have loaded the Amber module (e.g., 'module load amber') "
            "before running this script."
        )
    return amber_home

def check_amber_loaded():
    try:
        amber_path = get_amber_home_path()
    except EnvironmentError as e:
        print(f"Error: {e}")

        sys.exit(1)