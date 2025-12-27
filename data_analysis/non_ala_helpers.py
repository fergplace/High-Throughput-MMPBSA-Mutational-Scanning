
import pandas as pd
import numpy as np
from functools import  singledispatch
import matplotlib.pyplot as plt 
import seaborn as sns
from scipy import stats
from importlib import reload
import eval_funcs 
import plotting_funcs as pf 


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




def df_from_hkl(full_E_decomp : dict) -> pd.DataFrame:
    
    all_muts = []
    for mut_key in full_E_decomp.keys():
        mut_key_split = mut_key.split("_")
        mut =mut_key_split[0]
        all_muts.append(mut)
    
    all_muts= np.unique(np.array(all_muts, dtype=str))

    new_muts_dict = {key : {} for key in all_muts}

    for mut_key in full_E_decomp:
        mut_key_split = mut_key.split("_")
        mut =mut_key_split[0]
        radii = mut_key_split[1]
        model= mut_key_split[2]
        dg = float(full_E_decomp[mut_key]["DELTA TOTAL"][0])
        new_muts_dict[mut][(model + radii)]=  dg
        
    df = pd.DataFrame.from_dict(new_muts_dict, orient="index")
    return df



def df_from_hkl_std(full_E_decomp : dict) -> pd.DataFrame:
    
    all_muts = []
    for mut_key in full_E_decomp.keys():
        mut_key_split = mut_key.split("_")
        mut =mut_key_split[0]
        all_muts.append(mut)
    
    all_muts= np.unique(np.array(all_muts, dtype=str))

    new_muts_dict = {key : {} for key in all_muts}

    for mut_key in full_E_decomp:
        mut_key_split = mut_key.split("_")
        mut =mut_key_split[0]
        radii = mut_key_split[1]
        model= mut_key_split[2]
        dg = float(full_E_decomp[mut_key]["DELTA TOTAL"][0])
        new_muts_dict[mut][(model + radii + "dg")]=  dg
        new_muts_dict[mut][(model + radii + "std")]=  float(full_E_decomp[mut_key]["DELTA TOTAL"][1])
        
    df = pd.DataFrame.from_dict(new_muts_dict, orient="index")
    return df