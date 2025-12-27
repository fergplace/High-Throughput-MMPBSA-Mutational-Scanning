import pandas as pd
import numpy as np
from functools import  singledispatch
import matplotlib.pyplot as plt 
import seaborn as sns
from scipy import stats


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

@singledispatch
def relevant_exp_results( sims_done_arr, df_exp_data) : 
    tmp_dict= {}
    for key, value in df_exp_data["mutation"].to_dict().items() :
        if value[1:] in sims_done_arr :
            tmp_dict[key] = df_exp_data.loc[key]
            
    results_df = pd.DataFrame.from_dict(tmp_dict, orient="index")
    Rb_provided = results_df["bind_avg"].to_numpy()  
    R = 1.987*1e-3 #1e-3 for cal to kcal conversion 
    T =300  #room temp  
    K_D_app_wild = 3.9*1e-11 #from paper 

    K_D_app = 10**(Rb_provided)
    K_D_app_mut = K_D_app_wild /  K_D_app

    delta_g_deg_wild = -R * T * np.log(K_D_app_wild)
    del_del_g_mut = -R * T * np.log(K_D_app_mut)

    #flip sign to align with amber ALA scanning 
    del_del_g_total =  -(delta_g_deg_wild - del_del_g_mut) #-(del_del_g_mut- delta_g_deg_wild )

    #make it easier to read
    results_df.insert(6,"del_del_G", del_del_g_total )
    #results_df["del_del_G"] = del_del_g_total

    return results_df

@relevant_exp_results.register(pd.DataFrame)
def _( df_num_data, df_exp_data) : 
	sims_done_arr = df_num_data["mutation_idxed"].to_numpy()
	results_df = relevant_exp_results(sims_done_arr, df_exp_data)
	return results_df.sort_values(by = "site_SARS2")


def numerical_data_df_gen(csv_name :str):
    """ function for taking outputing df for numerical data given csv file name/path
    """
    df_num_data = [] 
    df_num_data = pd.read_csv(csv_name).drop(columns = "Unnamed: 0"	)
    mut_single_letter= []
    for idx, mut in enumerate(df_num_data["mutations"] ) : 
        mut_single_letter.append(three_letter_amino[mut ] )
    mut_single_letter = np.array(mut_single_letter)
    mutation_idxed = df_num_data["res_idx"].to_numpy(dtype=str) + mut_single_letter 
    df_num_data.insert( 2, "mutation_idxed" , mutation_idxed)
    df_num_data = df_num_data.sort_values(by="res_idx")

    return df_num_data


def get_data_helper(df_num_data,
                    del_g_wild_type_gb5 =-70.5334,
                    del_g_wild_type_gb66= -84.9824 #notew these are old values
                    ) :
    
    """ function used for getting the del del g values for the numerical data"""
    delta_g_gen_num_data_arr = df_num_data["delta_g_gen"].to_numpy()
    delta_g_gen_66_num_data_arr = df_num_data["delta_g_gen_66"].to_numpy()

    del_del_g_gen_num_data_arr  = delta_g_gen_num_data_arr- del_g_wild_type_gb5 
    
    del_del_g_gen_66_num_data_arr = delta_g_gen_66_num_data_arr - del_g_wild_type_gb66  
    
    return del_del_g_gen_num_data_arr,  del_del_g_gen_66_num_data_arr


def wildtype_mut_marker_helper(df_exp_data):
    """ funtion that returns all the wildtype markers e.g. 492L, used for checking the sim data to extract wildtype data (for our del del G calcs)
    """
    wildtype_df_tmp = df_exp_data[(df_exp_data["mutant"] == "*").to_numpy()]
    wildtype_mut_markers = wildtype_df_tmp["site_SARS2"].to_numpy(dtype=str) + wildtype_df_tmp["wildtype"].to_numpy() 
    return wildtype_mut_markers



def wildtype_df_gen(numerical_df, wild_markers):
    wildtype_marker= []
    for idx, row in numerical_df.iterrows():
        if row["mutation_idxed"] in wild_markers:
            wildtype_marker.append(True)
        else:
            wildtype_marker.append(False)
    wildtype_marker = np.array(wildtype_marker)
    numerical_df["wildtype_marker"] = wildtype_marker
    return numerical_df


def wildtype_calc(numerical_df):
    wildtype_rows = numerical_df[numerical_df["wildtype_marker"] == True]
    mean_gb5 = np.mean(wildtype_rows["delta_g_gen"].to_numpy())
    mean_gb66 = np.mean(wildtype_rows["delta_g_gen_66"].to_numpy())
    return {"mean_gb5": mean_gb5, "mean_gb66": mean_gb66}

def drop_wildtype(numerical_df):
    return numerical_df[numerical_df["wildtype_marker"] == False]


def add_del_del_g(numerical_df ) :
    
    "function that will return a numerical df with only the non wildtype rows and the del del g values added"
    wild_calc = wildtype_calc(numerical_df)
    num_df_no_wild = drop_wildtype(numerical_df).copy()
    
    # wild_calc["mean_gb5"] - num_df_no_wild["delta_g_gen"].to_numpy() 
    # wild_calc["mean_gb66"] - num_df_no_wild["delta_g_gen_66"].to_numpy()
    
    num_df_no_wild["del_del_g_gen"] = num_df_no_wild["delta_g_gen"].to_numpy() - wild_calc["mean_gb5"]
    num_df_no_wild["del_del_g_gen_66"] =  num_df_no_wild["delta_g_gen_66"].to_numpy() -wild_calc["mean_gb66"]
    return num_df_no_wild


def add_exp_data(df_num_data, df_exp_data) : 
    """ function that will add the experimental data to the numerical data"""
    del_del_g_exp  = relevant_exp_results(df_num_data, df_exp_data)["del_del_G"].to_numpy()
    df_num_data_updated = df_num_data.copy()
    df_num_data_updated["del_del_G_exp"] = del_del_g_exp
    return df_num_data_updated


def complete_num_df_gen(csv_path, df_exp_data):
    """ function that will take the numerical csv and return a df with del del G values added for both gb5 and gb66 and add the exp del del G values"""
    wild_markers = wildtype_mut_marker_helper(df_exp_data)
    df__num_data_raw = numerical_data_df_gen(csv_path)
    df_num_data_tmp = wildtype_df_gen(df__num_data_raw, wild_markers)
    df_num_data_del_del_g_added = add_del_del_g(df_num_data_tmp)
    df_num_data_exp_added = add_exp_data(df_num_data_del_del_g_added,df_exp_data )
    
    return df_num_data_exp_added


def mean_std_wildtype_df_gen(df_exp_data, csv_path):

    wild_markers = wildtype_mut_marker_helper(df_exp_data)
    df__num_data_raw = numerical_data_df_gen(csv_path)
    df_num_data_tmp = wildtype_df_gen(df__num_data_raw, wild_markers)
    df_num_data_tmp[df_num_data_tmp["wildtype_marker"] == True][["delta_g_gen","delta_g_gen_66"] ].to_numpy()
    
    wildtype_arr_mean = np.mean(df_num_data_tmp[df_num_data_tmp["wildtype_marker"] == True][["delta_g_gen","delta_g_gen_66"] ].to_numpy(), axis=0)
    wildtype_arr_std = np.std(df_num_data_tmp[df_num_data_tmp["wildtype_marker"] == True][["delta_g_gen","delta_g_gen_66"] ].to_numpy(), axis=0)
    # dumb, but works.... 
    wildtype_data = np.zeros(len(wildtype_arr_mean) + len(wildtype_arr_std))
    wildtype_data[0::2] = wildtype_arr_mean
    wildtype_data[1::2]= wildtype_arr_std
    wildtype_df = pd.DataFrame( columns=["gb5 mean","gb5 std",  "gb66 mean" ,"gb66 std"])
    wildtype_df.loc[0]=wildtype_data
    
    return wildtype_df