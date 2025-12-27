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

data_set_names = ["pb", "gb5_mbondi", "gb5_opt1", "gb5_opt5", "gb6_mbondi", "gb6_opt1", "gb6_opt5", "exp"]

colormap = plt.cm.get_cmap("viridis", len(data_set_names))
colors = {name: colormap(i) for i, name in enumerate(data_set_names)} #our color dict


def single_df_plotter(df : pd.DataFrame, radii : str, drop_gb5 = False, title_string = None) :
    num_rows = np.shape(df)[0]
    plt.figure(figsize=(10, 5))
    if drop_gb5 == False:
        plt.scatter(np.arange(0, num_rows) ,
                    df["del_del_g_gen"], 
                    label = f"gb5_{radii}", 
                    color = colors[f"gb5_{radii}"], edgecolors="black")
        
        correlation_coefficient, p_value = stats.pearsonr(df["del_del_g_gen"], df["del_del_G_exp"])
        
        mse = MSE(df["del_del_g_gen"], df["del_del_G_exp"])
        plt.plot([], [], ' ', label=f'MSE_5: {mse:.2f}')
        plt.plot([], [], ' ', label=f'corr coef _5: {correlation_coefficient:.2f}')
            
        
        
    plt.scatter(np.arange(0, num_rows) ,
                df["del_del_g_gen_66"], 
                label = f"gb66_{radii}", 
                color = colors[f"gb6_{radii}"], edgecolors="black")
    
    plt.scatter(np.arange(0, num_rows) ,
                df["del_del_G_exp"], 
                label = "exp_data", 
                color = colors[f"exp"], edgecolors="black")
    
    correlation_coefficient, p_value = stats.pearsonr(df["del_del_g_gen_66"], df["del_del_G_exp"])
    mse = MSE(df["del_del_g_gen_66"], df["del_del_G_exp"])
    plt.plot([], [], ' ', label=f'MSE_6: {mse:.2f}')
    plt.plot([], [], ' ', label=f'corr coef _6: {correlation_coefficient:.2f}')
    
    
    plt.xticks(np.arange(0, num_rows), df["mutation_idxed"].to_numpy(), rotation=45)
    plt.ylabel(r"$\Delta \Delta G$", fontsize=15)
    if title_string == None:
        plt.title(f"sim vs. exp for {radii}", fontsize=20)  
    else :
        plt.title(title_string, fontsize=20)
    plt.title(f"sim vs. exp for {radii}", fontsize=20)
    plt.legend(loc="upper right")
    plt.axhline(y= 0, color="grey", linestyle="--", zorder=0)
    
    plt.tight_layout()
    plt.savefig(f"single_piont_{radii}.png")
    
def MSE(arr1, arr2) :
    return np.mean((arr1 - arr2)**2)

def RMSE(arr1, arr2) :
    return np.sqrt(np.mean((arr1 - arr2)**2))

def std_error(arr1, arr2):
    return np.mean(arr1 - arr2)    #np.std((arr1 - arr2)) 


def mu(arr1, arr2):
    return np.mean(arr1 - arr2)

def multi_radii_plot(dfs :tuple, radii : tuple, drop_gb5= True, 
                     make_table = True):
    
    
    
    
    all_mse = [] # for the MSE and correlation coefficients, note we are only making the table for gb6 
    all_corr = []
    all_mean = []
    plt.figure(figsize=(10, 5))
    
    for df, radi in zip(dfs, radii):
        num_rows = np.shape(df)[0]
        if drop_gb5 == False and radi=="mbondi":
            plt.scatter(np.arange(0, num_rows) ,
                        df["del_del_g_gen"], 
                        label = f"gb5_{radi}", 
                        color = colors[f"gb5_{radi}"], edgecolors="black")
            
            correlation_coefficient, p_value = stats.pearsonr(df["del_del_g_gen"], df["del_del_G_exp"])
            mse = RMSE(df["del_del_g_gen"], df["del_del_G_exp"])
            mean = mu(df["del_del_g_gen"], df["del_del_G_exp"])
            all_mse.append(mse)
            all_mean.append(mean)
            all_corr.append(correlation_coefficient)
            plt.plot([], [], ' ', label=f'MSE_5: {mse:.2f}')
            plt.plot([], [], ' ', label=f'corr coef _5: {correlation_coefficient:.2f}')
                
            
            
        plt.scatter(np.arange(0, num_rows) ,
                    df["del_del_g_gen_66"], 
                    label = f"gb66_{radi}", 
                    color = colors[f"gb6_{radi}"], edgecolors="black")
        
        
        
        correlation_coefficient, p_value = stats.pearsonr(df["del_del_g_gen_66"], df["del_del_G_exp"])
        mse = RMSE(df["del_del_g_gen_66"], df["del_del_G_exp"])
        mean = mu(df["del_del_g_gen_66"], df["del_del_G_exp"])
        
        all_mse.append(mse)
        all_mean.append(mean)
        all_corr.append(correlation_coefficient)
        
        # plt.plot([], [], ' ', label=f'MSE_6: {mse:.2f}')
        # plt.plot([], [], ' ', label=f'corr coef _6: {correlation_coefficient:.2f}')
        
        
    
    plt.scatter(np.arange(0, num_rows) ,
                    df["del_del_G_exp"], 
                    label = "exp_data", 
                    color = colors[f"exp"], edgecolors="black")
    
    plt.xticks(np.arange(0, num_rows), df["mutation_idxed"].to_numpy(), rotation=45)
    plt.ylabel(r"$\Delta \Delta G$", fontsize=15)
    plt.title(f"sim vs. exp for different radii", fontsize=20)
    plt.legend(loc="upper right")
    plt.axhline(y= 0, color="grey", linestyle="--", zorder=0)
    plt.tight_layout()
    plt.savefig("single_point_all_radii.png")
    radii = ["gb5 mbondi","gb66 mbondi", "opt1", "opt5"]
    if make_table == True:
        table = pd.DataFrame({"Radii": radii,
                              "mean": all_mean, 
                              "RMSE": all_mse, 
                              "Correlation Coefficient": all_corr})
        return table
    
    
    