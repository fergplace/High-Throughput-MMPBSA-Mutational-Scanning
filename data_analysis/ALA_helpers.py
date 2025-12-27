import pandas as pd
import numpy as np
from functools import  singledispatch
import matplotlib.pyplot as plt 
import seaborn as sns
from scipy import stats
import eval_funcs 
import plotting_funcs as pf
from scipy.stats import skew
import numpy.typing as npt
import copy 

"""general code for plottig and table gen. This code here is a bit funcky, but it mainly serves to populate tables and figures with correct labels
"""
def df_mapper_with_pb(df_all, radii):
	df = df_all[df_all["radii"] == radii].copy()
	new_mut_names = [] 
	mutation_idx_key_helper = [] 
	for mut in df["mutations"].to_numpy() :
		new_mut_names.append(mut)
		mutation_idx_key_helper.append(mut[1:])
	df["mutations"] = np.array(new_mut_names)
	df.insert(1, "mutation_idxed",  np.array(mutation_idx_key_helper )) 
	df.drop(columns=["delta_delta_g_pb_error", "delta_delta_g_gen_66_error" , "delta_delta_g_gen_error", "radii"], inplace=True)
	df = df.rename(columns={"delta_delta_g_pb" : f"ddg_pb_{radii}", 
						"delta_delta_g_gen_66" : f"ddg_gb_66_{radii}",
					"delta_delta_g_gen" : f"ddg_gb_5_{radii}"})
	df = df.sort_values(by="mutation_idxed", ignore_index=True)
	return df 

def sim_data_simple_df_gen_with_pb(path : str, radii ):
    df_raw = pd.read_csv(path, index_col=0)
    df = df_raw.copy() #just make a deep copy 
    new_mut_names = [] 
    mutation_idx_key_helper = [] 
    for mut in df["mutations"].to_numpy() :
        tmp_str = mut.split("_")[0]
        new_mut_names.append(tmp_str)
        mutation_idx_key_helper.append(tmp_str[1:])

    df["mutations"] = np.array(new_mut_names)
    df.insert(1, "mutation_idxed",  np.array(mutation_idx_key_helper )) 

    # now drop the errors for now, we can change this later, but for now we are just interested in the del del 

    df.drop(columns=["delta_delta_g_pb_error", "delta_delta_g_gen_66_error" , "delta_delta_g_gen_error"], inplace=True)
    df = df.rename(columns={"delta_delta_g_pb" : f"ddg_pb_{radii}", 
                        "delta_delta_g_gen_66" : f"ddg_gb_66_{radii}",
					"delta_delta_g_gen" : f"ddg_gb_5_{radii}"})
    df = df.sort_values(by="mutation_idxed", ignore_index=True)
    return df

def sim_data_simple_df_gen_without_pb(path : str, radii ):
    df_raw = pd.read_csv(path)
    df = df_raw.copy() #just make a deep copy 
    new_mut_names = [] 
    mutation_idx_key_helper = [] 
    for mut in df["mutations"].to_numpy() :
        tmp_str = mut.split("_")[0]
        new_mut_names.append(tmp_str)
        mutation_idx_key_helper.append(tmp_str[1:])

    df["mutations"] = np.array(new_mut_names)
    df.insert(1, "mutation_idxed",  np.array(mutation_idx_key_helper )) 

    # now drop the errors for now, we can change this later, but for now we are just interested in the del del 
    cols_to_drop = _df_drop_helper(df)


    df.drop(columns=cols_to_drop, inplace=True)
    
    df = df.sort_values(by="mutation_idxed", ignore_index=True)
    return df

def _condense_col_names(column_names :  npt.NDArray[np.str_]):
    new_names = [] 
    for col in column_names:
        col_copy = copy.copy(col)
        if col_copy.startswith("delta_delta_g_gen_66") :
                new_names.append(col_copy.replace("delta_delta_g_gen", "ddg_gb"))
        elif col_copy.startswith("delta_delta_g_gen") :
                new_names.append(col_copy.replace("delta_delta_g_gen", "ddg_gb_5"))
        elif col_copy.startswith("delta_delta_g_pb") :
                new_names.append(col_copy.replace("delta_delta_g_pb", "ddg_pb"))
        else: 
                new_names.append(col_copy) 
    return new_names

def _default_mut_str_splitter(mut_str : str) -> str:
    return mut_str.split("_")[0]


def _add_mutation_idxed(df :pd.DataFrame, 
                        default_mut_str_splitter = _default_mut_str_splitter):
       df_copy= df.copy() #make a copy of the df aimsing for immuatble states
       new_mut_names = [] 
       mutation_idx_key_helper = [] 
       for mut in df["mutations"].to_numpy() :
              tmp_str = _default_mut_str_splitter(str(mut))
              new_mut_names.append(tmp_str)
              mutation_idx_key_helper.append(tmp_str[1:])

       df_copy["mutations"] = np.array(new_mut_names)
       df_copy.insert(1, "mutation_idxed",  np.array(mutation_idx_key_helper )) 
       return df_copy 


def _df_drop_helper(df : pd.DataFrame, str_to_drop_by : str = "_error"):
       cols = df.columns.to_numpy()
       cols_to_drop = [col for col in cols if f"{str_to_drop_by}" in col]
       return cols_to_drop



def check_valid_dfs(df1, df2):
    "simple way to check if the dfs are ordered in the same way"
    return np.all(df1["mutation_idxed"].to_numpy() == df2["mutation_idxed"].to_numpy() )

def full_df_output_generator(dfs : tuple):
    """ function that takes in a tuple of dfs and a tuple of radii and outputs a df with all the data. Note first df should contain the pb data (mbondi will do)
    """
    df = dfs[0].copy() #copy first df, this is the one we will be modifying 
    
    for idx, tmp_df in enumerate(dfs[1:]):
        if not check_valid_dfs(dfs[0], tmp_df):
            raise ValueError("dfs are not ordered in the same way.. fix it")
        else:
            for col in tmp_df.columns.to_numpy() : 
                if "_".join(col.split("_")[:-1]) in ["ddg_gb_5", "ddg_gb_66"]:
                    df[col] = tmp_df[col]
                    
    return df


def add_exp_data(num_df, exp_df ) :
    relevant_exp_df = eval_funcs.relevant_exp_results( num_df, exp_df)
    if not np.all(relevant_exp_df["mutation"].to_numpy() == num_df["mutations"].to_numpy() ) :
        raise ValueError("exp data and numerical data are not ordered in the same way")
    
    num_df["ddg_exp"] = relevant_exp_df["del_del_G"].to_numpy()
    
    return num_df

#array for binding site
binding_site =np.array( ['K417A', 'N439A', 'Y449A', 'Y453A', 'L455A', 'F456A',
	'F486A', 'N487A', 'Y489A', 'Q493A', 'S494A', 'Q498A',
	'T500A', 'N501A', 'V503A', 'Y505A'])

def add_binding_site_bool(df) :
    binding_bool_arr = np.array([(site in binding_site) for site in df["mutations"]] , dtype=int )
    df_bd_site_added = df.copy()
    df_bd_site_added["bd_site"] = binding_bool_arr
    return df_bd_site_added


data_set_names = ["pb", "gb5_mbondi2", "gb5_opt1", "gb5_opt5", "gb6_mbondi2", "gb6_opt1", "gb6_opt5", "exp"]

colormap = plt.cm.get_cmap("viridis", len(data_set_names))
colors = {name: colormap(i) for i, name in enumerate(data_set_names)} #our color dict


#plotting funcs

def full_df_ala_plotter(df, full=True, 
                        pb=False,
                        gb_mbondi = False,
                        save_fig=False,
                        fig_title : str = None):
    
    num_rows = np.shape(df)[0]
    scatter_x_arr = np.arange(0, num_rows)
    plt.figure(figsize=(10, 5))
    plt.scatter(scatter_x_arr, df["ddg_gb_66_mbondi2"].to_numpy(), color=colors["gb6_mbondi2"], edgecolors="black", label="gb6_mbondi2")

    plt.scatter(scatter_x_arr, df["ddg_gb_66_opt1"], color=colors["gb6_opt1"], edgecolors="black", label="gb6_opt1")


    plt.scatter(scatter_x_arr, df["ddg_gb_66_opt5d"], color=colors["gb6_opt5"], edgecolors="black", label="gb6_opt5")


    plt.scatter(scatter_x_arr, df["ddg_exp"], color=colors["exp"], edgecolors="black", label="exp")

    
    if full == True:
        tick_positions = np.arange(0, len(df["mutations"]), 5)
        tick_labels = [df["mutations"].iloc[pos] for pos in tick_positions]
        tick_labels = [df["mutations"].iloc[pos][1:-1] for pos in tick_positions]
        plt.xticks(tick_positions, tick_labels,rotation=45,  fontsize = 8)
        plt.xlabel("Residue Number", fontsize = 15)
    else :
        tick_positions = np.arange(0, len(df["mutations"]), 1)
        tick_labels = [df["mutations"].iloc[pos] for pos in tick_positions]
        tick_labels = [df["mutations"].iloc[pos][1:-1] for pos in tick_positions]
        plt.xticks(tick_positions, tick_labels, rotation=45)
        plt.xlabel("Residue",  fontsize = 20)
    if pb == True:
        plt.scatter(scatter_x_arr, df["ddg_pb_mbondi2"], color=colors["pb"], edgecolors="black", label="pb")

    if gb_mbondi == True :
        plt.scatter(scatter_x_arr, df["ddg_gb_5_mbondi2"], color=colors["gb5_mbondi2"] , label="gb5_mbondi2", marker="x") # edgecolors="black"


    plt.yticks(fontsize= 15)


    plt.legend()
    plt.title("All ALA mutations over different radii for gb6")
    plt.ylabel(r"$\Delta \Delta G$")
    
    if save_fig ==True:
        if fig_title is not None:
            plt.savefig(f"{fig_title}.png")
        else:
            raise ValueError("need fig title")
        
        

def df_fist_pass_plotter(df, full=True, 
                        pb=False,
                        gb_mbondi = False,
                        save_fig=False,
                        fig_title : str = None):
    
    num_rows = np.shape(df)[0]
    scatter_x_arr = np.arange(0, num_rows)
    plt.figure(figsize=(10, 5))
    plt.scatter(scatter_x_arr, df["ddg_gb_66_mbondi2"].to_numpy(), color=colors["gb6_mbondi2"], edgecolors="black", label="gb6_mbondi2")

    plt.scatter(scatter_x_arr, df["ddg_exp"], color=colors["exp"], edgecolors="black", label="exp")

    
    if full == True:
        tick_positions = np.arange(0, len(df["mutations"]), 5)
        tick_labels = [df["mutations"].iloc[pos] for pos in tick_positions]
        tick_labels = [df["mutations"].iloc[pos][1:-1] for pos in tick_positions]
        plt.xticks(tick_positions, tick_labels,rotation=45,  fontsize = 8)
        plt.xlabel("Residue Number", fontsize = 15)
    else :
        tick_positions = np.arange(0, len(df["mutations"]), 1)
        tick_labels = [df["mutations"].iloc[pos] for pos in tick_positions]
        tick_labels = [df["mutations"].iloc[pos][1:-1] for pos in tick_positions]
        plt.xticks(tick_positions, tick_labels, rotation=45)
        plt.xlabel("Residue",  fontsize = 20)
    if pb == True:
        plt.scatter(scatter_x_arr, df["ddg_pb_mbondi2"], color=colors["pb"], edgecolors="black", label="pb")

    if gb_mbondi == True :
        plt.scatter(scatter_x_arr, df["ddg_gb_5_mbondi2"], color=colors["gb5_mbondi2"] , label="gb5_mbondi2", marker="x") # edgecolors="black"


    plt.yticks(fontsize= 15)

    fig_title_as_label = " ".join( fig_title.split("_") )
    plt.legend()
    plt.title(f" {  fig_title_as_label}")
    plt.ylabel(r"$\Delta \Delta G$")
    
    if save_fig ==True:
        if fig_title is not None:
            plt.savefig(f"{fig_title}.png")
        else:
            raise ValueError("need fig title")
        
        
def first_pass_table_gen(df) :
    exp_vals = df["ddg_exp"].to_numpy()
    mse_exp = np.array( [ pf.RMSE(exp_vals, df["ddg_gb_5_mbondi2"]) , 
                         pf.RMSE(exp_vals, df["ddg_gb_66_mbondi2"]),
                         pf.RMSE(exp_vals, df["ddg_pb_mbondi2"]) ] )
    
    
    mse_pb  = np.array( [ pf.RMSE(df["ddg_pb_mbondi2"], df["ddg_gb_5_mbondi2"]) , 
                         pf.RMSE(df["ddg_pb_mbondi2"], df["ddg_gb_66_mbondi2"]) , 0.0 ] )

    
    df_binding = df[df["bd_site"] == True ].copy()
 
    mse_exp_binding = np.array( 
                    [ 
                    pf.RMSE(df_binding["ddg_exp"], 
                           df_binding["ddg_gb_5_mbondi2"]) , 
                    pf.RMSE(df_binding["ddg_exp"], 
                           df_binding["ddg_gb_66_mbondi2"]),
                    pf.RMSE(df_binding["ddg_exp"], 
                           df_binding["ddg_pb_mbondi2"]) 
                    ] )


    mse_pb_binding  = np.array( 
                    [ 
                pf.RMSE(df_binding["ddg_pb_mbondi2"],
                       df_binding["ddg_gb_5_mbondi2"]) , 
                pf.RMSE(df_binding["ddg_pb_mbondi2"], 
                       df_binding["ddg_gb_66_mbondi2"]) ,
                    0.0 
                    ] )

    #non binding 
    df_non_binding = df[df["bd_site"] == False ].copy()
    mse_exp_non_binding = np.array( 
                    [ 
                    pf.RMSE(df_non_binding["ddg_exp"], 
                           df_non_binding["ddg_gb_5_mbondi2"]) , 
                    pf.RMSE(df_non_binding["ddg_exp"], 
                           df_non_binding["ddg_gb_66_mbondi2"]),
                    pf.RMSE(df_non_binding["ddg_exp"], 
                           df_non_binding["ddg_pb_mbondi2"]) 
                    ] )



    mse_pb_non_binding  = np.array( 
                    [ 
                pf.RMSE(df_non_binding["ddg_pb_mbondi2"],
                       df_non_binding["ddg_gb_5_mbondi2"]) , 
                pf.RMSE(df_non_binding["ddg_pb_mbondi2"], 
                       df_non_binding["ddg_gb_66_mbondi2"]) ,
                    0.0 
                    ] )
    
    ##### 

    tmp_df = pd.DataFrame.from_dict( 
                            { 
                                "Model" : ["gb5 mbondi2", "gb 66 mbondi2", "pb"] ,
                                "RMSE exp" : mse_exp , 
                                "RMSE exp N-BS" : mse_exp_non_binding,
                                "RMSE exp binding" : mse_exp_binding,
                                "RMSE Pb" : mse_pb ,
                                "RMSE pb N-BS" : mse_pb_non_binding,
                                "RMSE pb binding" : mse_pb_binding
                            } )
    
        
    tmp_df.set_index('Model', inplace=True)
    return tmp_df


def corr_plot_first_pass(full_data_df_bd_site, save_fig = False):
    f, ax = plt.subplots(2, 3, figsize=(12,6), height_ratios=[0.9, 0.1])

    corr = full_data_df_bd_site[["ddg_exp" , "ddg_pb_mbondi2", "ddg_gb_5_mbondi2", "ddg_gb_66_mbondi2"]].corr()
    corr_plot_labels = ["PB", "GB5 mbondi2", "GB 66 mbondi2"]
    sns.heatmap(corr, 
                xticklabels=[0, 1, 2, 3],
                yticklabels=[0, 1, 2, 3], cmap='coolwarm', cbar=False, ax = ax[0,0])
    ax[0, 0].set_title("Full RBD")

    corr = full_data_df_bd_site[["ddg_exp" , "ddg_pb_mbondi2", "ddg_gb_5_mbondi2", "ddg_gb_66_mbondi2"]][full_data_df_bd_site["bd_site"] == True ].corr()
    sns.heatmap(corr, 
                xticklabels=[0, 1, 2, 3], yticklabels=False, cmap='coolwarm', cbar=False,  ax = ax[0,1])
    ax[0, 1].set_title("Binding Site")

    corr = full_data_df_bd_site[["ddg_exp" ,"ddg_pb_mbondi2", "ddg_gb_5_mbondi2", "ddg_gb_66_mbondi2"]][full_data_df_bd_site["bd_site"] == False ].corr()
    gs = ax[1, 2].get_gridspec()
    axbig = f.add_subplot(gs[1,:2])
    sns.heatmap(corr, 
                xticklabels=[0, 1, 2, 3], yticklabels=False,cmap='coolwarm', cbar=True,  ax = ax[0,2], cbar_ax=axbig, cbar_kws={"orientation" : "horizontal", "label" : "Corelation Coef"})
    ax[0, 2].set_title("Non-Binding Site")
    for ax in ax[1, :3]:
        ax.remove()
    legend_ax = f.add_subplot(gs[1, 2])
    legend_ax.axis("off") 
    legend_labels = [
        "0: exp",
        "1: PB",
        "2: GB 5 mbondi2",
        "3: GB 66 mbondi2",
    ]
    legend_handles = [plt.Line2D([0], [0], color='none', label=label) for label in legend_labels]
    legend_ax.legend(handles=legend_handles, loc="upper center", fontsize=12, frameon=True, edgecolor="black")
    f.suptitle("Correlation Plots for mbondi2", fontsize=20)
    f.tight_layout()

    if save_fig == True:
        f.savefig("Correlation_Plots_for_mbondi2")
    return corr
        
        
def radii_table_gen(df) :
    exp_vals = df["ddg_exp"].to_numpy()
    mse_exp = np.array( [ pf.MSE(exp_vals, df["ddg_gb_5_mbondi2"]) , 
                         pf.MSE(exp_vals, df["ddg_gb_66_mbondi2"]),
                         pf.MSE(exp_vals, df["ddg_gb_66_opt1"]),
                         pf.MSE(exp_vals, df["ddg_gb_66_opt5d"]),
                         pf.MSE(exp_vals, df["ddg_pb_mbondi2"]) ] )
    
    
    mse_pb  = np.array( [ 
                         pf.MSE(df["ddg_pb_mbondi2"], df["ddg_gb_5_mbondi2"]) , 
                         pf.MSE(df["ddg_pb_mbondi2"], df["ddg_gb_66_mbondi2"]) , 
                         pf.MSE(df["ddg_pb_mbondi2"], df["ddg_gb_66_opt1"]) , 
                         pf.MSE(df["ddg_pb_mbondi2"], df["ddg_gb_66_opt5d"]) , 
                         0.0 ] )

    
    df_binding = df[df["bd_site"] == True ].copy()
    mse_exp_binding = np.array( 
                    [ 
                    pf.MSE(df_binding["ddg_exp"], 
                           df_binding["ddg_gb_5_mbondi2"]) , 
                    
                    pf.MSE(df_binding["ddg_exp"], 
                           df_binding["ddg_gb_66_mbondi2"]),
                    pf.MSE(df_binding["ddg_exp"], 
                           df_binding["ddg_gb_66_opt1"]),
                    pf.MSE(df_binding["ddg_exp"], 
                           df_binding["ddg_gb_66_opt5d"]) ,
                    
                    pf.MSE(df_binding["ddg_exp"], 
                           df_binding["ddg_pb_mbondi2"]) 
                    
                    ] )


    mse_pb_binding  = np.array( 
                    [ 
                pf.MSE(df_binding["ddg_pb_mbondi2"],
                       df_binding["ddg_gb_5_mbondi2"]) , 
                
                pf.MSE(df_binding["ddg_pb_mbondi2"], 
                       df_binding["ddg_gb_66_mbondi2"]) ,
                
                pf.MSE(df_binding["ddg_pb_mbondi2"], 
                           df_binding["ddg_gb_66_opt1"]),
                pf.MSE(df_binding["ddg_pb_mbondi2"], 
                           df_binding["ddg_gb_66_opt5d"]) ,
                    0.0 
                    ] )

    #non binding 
    df_non_binding = df[df["bd_site"] == False ].copy()
    mse_exp_non_binding = np.array( 
                    [ 
                    pf.MSE(df_non_binding["ddg_exp"], 
                           df_non_binding["ddg_gb_5_mbondi2"]) , 
                    pf.MSE(df_non_binding["ddg_exp"], 
                           df_non_binding["ddg_gb_66_mbondi2"]),
                    pf.MSE(df_non_binding["ddg_exp"], 
                           df_non_binding["ddg_gb_66_opt1"]),
                    pf.MSE(df_non_binding["ddg_exp"], 
                           df_non_binding["ddg_gb_66_opt5d"]) ,
                    
                    pf.MSE(df_non_binding["ddg_exp"], 
                           df_non_binding["ddg_pb_mbondi2"]) 
                    
                    ] )



    mse_pb_non_binding  = np.array( 
                    [ 
                pf.MSE(df_non_binding["ddg_pb_mbondi2"],
                       df_non_binding["ddg_gb_5_mbondi2"]) ,
                 
                pf.MSE(df_non_binding["ddg_pb_mbondi2"], 
                       df_non_binding["ddg_gb_66_mbondi2"]) ,
                
                pf.MSE(df_non_binding["ddg_pb_mbondi2"], 
                        df_non_binding["ddg_gb_66_opt1"]),
                
                pf.MSE(df_non_binding["ddg_pb_mbondi2"], 
                        df_non_binding["ddg_gb_66_opt5d"]) ,
                    
                    0.0 
                    ] )
    
    ##### 

    tmp_df = pd.DataFrame.from_dict( 
                            { 
                                "Model" : ["gb5 mbondi2", 
                                           "gb 66 mbondi2",
                                           "gb 66 op1" ,
                                           "gb 66 op5d" ,
                                           "pb"] ,
                                "MSE exp" : mse_exp , 
                                "MSE exp N-BS" : mse_exp_non_binding,
                                "MSE exp binding" : mse_exp_binding,
                                "MSE Pb" : mse_pb ,
                                "MSE pb N-BS" : mse_pb_non_binding,
                                "MSE pb binding" : mse_pb_binding
                            } )
    
        
    tmp_df.set_index('Model', inplace=True)
    return tmp_df




def corr_plot_radii(full_data_df_bd_site, save_fig = False):
    f, ax = plt.subplots(2, 3, figsize=(12,5), height_ratios=[0.9, 0.1])

    desired_cols = ["ddg_exp" , 
                    "ddg_pb_mbondi2", 
                    "ddg_gb_5_mbondi2", 
                    "ddg_gb_66_mbondi2", 
                    "ddg_gb_66_opt1", 
                    "ddg_gb_66_opt5d"]
    corr = full_data_df_bd_site[desired_cols].corr()
    labels = np.arange(0, len(desired_cols))
    sns.heatmap(corr, 
                xticklabels=labels,
                yticklabels=labels, cmap='coolwarm', cbar=False, ax = ax[0,0], square=True)
    ax[0, 0].set_title("Full RBD")

    corr = full_data_df_bd_site[desired_cols][full_data_df_bd_site["bd_site"] == True ].corr()
    
    
    sns.heatmap(corr, 
                xticklabels=labels, yticklabels=False, cmap='coolwarm', cbar=False,  ax = ax[0,1], square=True)
    ax[0, 1].set_title("Binding Site")

    corr = full_data_df_bd_site[desired_cols][full_data_df_bd_site["bd_site"] == False ].corr()
    gs = ax[1, 2].get_gridspec()
    axbig = f.add_subplot(gs[1,:2])
    sns.heatmap(corr, 
                xticklabels=labels, yticklabels=False,cmap='coolwarm', cbar=True,  ax = ax[0,2], cbar_ax=axbig, cbar_kws={"orientation" : "horizontal", "label" : "Corelation Coef"}, square=True)
    ax[0, 2].set_title("Non-Binding Site")
    for ax in ax[1, :3]:
        ax.remove()
    legend_ax = f.add_subplot(gs[1, 2])
    legend_ax.axis("off") 
    legend_labels = [
        "0: exp",
        "1: PB",
        "2: GB 5 mbondi2",
        "3: GB 66 mbondi2",
        "4: GB 66 opt1",
        "5: GB 66 opt5d"
    ]
    legend_handles = [plt.Line2D([0], [0], color='none', label=label) for label in legend_labels]
    legend_ax.legend(handles=legend_handles, loc="upper center", fontsize=12, frameon=True, edgecolor="black", ncol=2)
    f.suptitle("Correlation Plots for mbondi2", fontsize=20)
    f.tight_layout()

    if save_fig == True:
        f.savefig("Correlation_Plots_for_all_radii")
    
    return corr
    


def rmse_std_table_gen(df) :
    exp_vals = df["ddg_exp"].to_numpy()
    mse_exp = np.array( [ pf.RMSE(exp_vals, df["ddg_gb_5_mbondi2"]) , 
                         pf.RMSE(exp_vals, df["ddg_gb_66_mbondi2"]),
                         pf.RMSE(exp_vals, df["ddg_pb_mbondi2"]) ] )
    std_exp = np.array( [pf.mu(exp_vals, df["ddg_gb_5_mbondi2"]) , 
                        pf.mu(exp_vals, df["ddg_gb_66_mbondi2"]),
                        pf.mu(exp_vals, df["ddg_pb_mbondi2"]) ] )
    
    mse_pb  = np.array( [ pf.RMSE(df["ddg_pb_mbondi2"], df["ddg_gb_5_mbondi2"]) , 
                         pf.RMSE(df["ddg_pb_mbondi2"], df["ddg_gb_66_mbondi2"]) , 0.0 ] )
    std_pb  = np.array( [pf.mu(df["ddg_pb_mbondi2"], df["ddg_gb_5_mbondi2"]) , 
                           pf.mu(df["ddg_pb_mbondi2"], df["ddg_gb_66_mbondi2"]) , 0.0 ] )
    
    df_binding = df[df["bd_site"] == True ].copy()

    mse_exp_binding = np.array( 
                    [ 
                    pf.RMSE(df_binding["ddg_exp"], 
                           df_binding["ddg_gb_5_mbondi2"]) , 
                    pf.RMSE(df_binding["ddg_exp"], 
                           df_binding["ddg_gb_66_mbondi2"]),
                    pf.RMSE(df_binding["ddg_exp"], 
                           df_binding["ddg_pb_mbondi2"]) 
                    ] )

    std_exp_binding = np.array( 
                [ 
               pf.mu(df_binding["ddg_exp"], 
                        df_binding["ddg_gb_5_mbondi2"]) , 
               pf.mu(df_binding["ddg_exp"], 
                        df_binding["ddg_gb_66_mbondi2"]),
               pf.mu(df_binding["ddg_exp"], 
                        df_binding["ddg_pb_mbondi2"]) 
                ] )
    
    mse_pb_binding  = np.array( 
                    [ 
                pf.RMSE(df_binding["ddg_pb_mbondi2"],
                       df_binding["ddg_gb_5_mbondi2"]) , 
                pf.RMSE(df_binding["ddg_pb_mbondi2"], 
                       df_binding["ddg_gb_66_mbondi2"]) ,
                    0.0 
                    ] )
    
    std_pb_binding  = np.array( 
                [ 
           pf.mu(df_binding["ddg_pb_mbondi2"],
                    df_binding["ddg_gb_5_mbondi2"]) , 
           pf.mu(df_binding["ddg_pb_mbondi2"], 
                    df_binding["ddg_gb_66_mbondi2"]) ,
                0.0 
                ] )

    #non binding 
    df_non_binding = df[df["bd_site"] == False ].copy()
    mse_exp_non_binding = np.array( 
                    [ 
                    pf.RMSE(df_non_binding["ddg_exp"], 
                           df_non_binding["ddg_gb_5_mbondi2"]) , 
                    pf.RMSE(df_non_binding["ddg_exp"], 
                           df_non_binding["ddg_gb_66_mbondi2"]),
                    pf.RMSE(df_non_binding["ddg_exp"], 
                           df_non_binding["ddg_pb_mbondi2"]) 
                    ] )
    std_exp_non_binding = np.array( 
                    [ 
                   pf.mu(df_non_binding["ddg_exp"], 
                            df_non_binding["ddg_gb_5_mbondi2"]) , 
                   pf.mu(df_non_binding["ddg_exp"], 
                            df_non_binding["ddg_gb_66_mbondi2"]),
                   pf.mu(df_non_binding["ddg_exp"], 
                            df_non_binding["ddg_pb_mbondi2"]) 
                    ] )


    mse_pb_non_binding  = np.array( 
                    [ 
                pf.RMSE(df_non_binding["ddg_pb_mbondi2"],
                       df_non_binding["ddg_gb_5_mbondi2"]) , 
                pf.RMSE(df_non_binding["ddg_pb_mbondi2"], 
                       df_non_binding["ddg_gb_66_mbondi2"]) ,
                    0.0 
                    ] )
    std_pb_non_binding  = np.array( 
                [ 
           pf.mu(df_non_binding["ddg_pb_mbondi2"],
                    df_non_binding["ddg_gb_5_mbondi2"]) , 
           pf.mu(df_non_binding["ddg_pb_mbondi2"], 
                    df_non_binding["ddg_gb_66_mbondi2"]) ,
                0.0 
                ] )

    ##### 

    tmp_df = pd.DataFrame.from_dict( 
                            { 
                                "Model" : ["gb5 mbondi2", "gb 66 mbondi2", "pb"] ,
                                "RMSE exp" : mse_exp , 
                                "std exp" : std_exp, 
                                
                                "RMSE exp N-BS" : mse_exp_non_binding,
                                "std exp N-BS" : std_exp_non_binding,
                                
                                "RMSE exp binding" : mse_exp_binding,
                                "std exp binding" : std_exp_binding,
                                
                                "RMSE Pb" : mse_pb ,
                                "std Pb" : std_pb ,
                                
                                "RMSE pb N-BS" : mse_pb_non_binding,
                                "std pb N-BS" : std_pb_non_binding,
                                
                                "RMSE pb binding" : mse_pb_binding,
                                "std pb binding" : std_pb_binding
                            } )
    
        
    tmp_df.set_index('Model', inplace=True)
    return tmp_df



def mu_rmse_table_gen(df) :
    exp_vals = df["ddg_exp"].to_numpy()
    mse_exp = np.array( [ pf.RMSE(exp_vals, df["ddg_gb_5_mbondi2"]) , 
                         pf.RMSE(exp_vals, df["ddg_gb_66_mbondi2"]),
                         pf.RMSE(exp_vals, df["ddg_pb_mbondi2"]) ] )
    std_exp = np.array( [pf.mu(exp_vals, df["ddg_gb_5_mbondi2"]) , 
                        pf.mu(exp_vals, df["ddg_gb_66_mbondi2"]),
                        pf.mu(exp_vals, df["ddg_pb_mbondi2"]) ] )
    
    mse_pb  = np.array( [ pf.RMSE(df["ddg_pb_mbondi2"], df["ddg_gb_5_mbondi2"]) , 
                         pf.RMSE(df["ddg_pb_mbondi2"], df["ddg_gb_66_mbondi2"]) , 0.0 ] )
    std_pb  = np.array( [pf.mu(df["ddg_pb_mbondi2"], df["ddg_gb_5_mbondi2"]) , 
                           pf.mu(df["ddg_pb_mbondi2"], df["ddg_gb_66_mbondi2"]) , 0.0 ] )
    
    df_binding = df[df["bd_site"] == True ].copy()
    mse_exp_binding = np.array( 
                    [ 
                    pf.RMSE(df_binding["ddg_exp"], 
                           df_binding["ddg_gb_5_mbondi2"]) , 
                    pf.RMSE(df_binding["ddg_exp"], 
                           df_binding["ddg_gb_66_mbondi2"]),
                    pf.RMSE(df_binding["ddg_exp"], 
                           df_binding["ddg_pb_mbondi2"]) 
                    ] )

    std_exp_binding = np.array( 
                [ 
               pf.mu(df_binding["ddg_exp"], 
                        df_binding["ddg_gb_5_mbondi2"]) , 
               pf.mu(df_binding["ddg_exp"], 
                        df_binding["ddg_gb_66_mbondi2"]),
               pf.mu(df_binding["ddg_exp"], 
                        df_binding["ddg_pb_mbondi2"]) 
                ] )
    
    mse_pb_binding  = np.array( 
                    [ 
                pf.RMSE(df_binding["ddg_pb_mbondi2"],
                       df_binding["ddg_gb_5_mbondi2"]) , 
                pf.RMSE(df_binding["ddg_pb_mbondi2"], 
                       df_binding["ddg_gb_66_mbondi2"]) ,
                    0.0 
                    ] )
    
    std_pb_binding  = np.array( 
                [ 
           pf.mu(df_binding["ddg_pb_mbondi2"],
                    df_binding["ddg_gb_5_mbondi2"]) , 
           pf.mu(df_binding["ddg_pb_mbondi2"], 
                    df_binding["ddg_gb_66_mbondi2"]) ,
                0.0 
                ] )

    #non binding 
    df_non_binding = df[df["bd_site"] == False ].copy()
    mse_exp_non_binding = np.array( 
                    [ 
                    pf.RMSE(df_non_binding["ddg_exp"], 
                           df_non_binding["ddg_gb_5_mbondi2"]) , 
                    pf.RMSE(df_non_binding["ddg_exp"], 
                           df_non_binding["ddg_gb_66_mbondi2"]),
                    pf.RMSE(df_non_binding["ddg_exp"], 
                           df_non_binding["ddg_pb_mbondi2"]) 
                    ] )
    std_exp_non_binding = np.array( 
                    [ 
                   pf.mu(df_non_binding["ddg_exp"], 
                            df_non_binding["ddg_gb_5_mbondi2"]) , 
                   pf.mu(df_non_binding["ddg_exp"], 
                            df_non_binding["ddg_gb_66_mbondi2"]),
                   pf.mu(df_non_binding["ddg_exp"], 
                            df_non_binding["ddg_pb_mbondi2"]) 
                    ] )


    mse_pb_non_binding  = np.array( 
                    [ 
                pf.RMSE(df_non_binding["ddg_pb_mbondi2"],
                       df_non_binding["ddg_gb_5_mbondi2"]) , 
                pf.RMSE(df_non_binding["ddg_pb_mbondi2"], 
                       df_non_binding["ddg_gb_66_mbondi2"]) ,
                    0.0 
                    ] )
    std_pb_non_binding  = np.array( 
                [ 
           pf.mu(df_non_binding["ddg_pb_mbondi2"],
                    df_non_binding["ddg_gb_5_mbondi2"]) , 
           pf.mu(df_non_binding["ddg_pb_mbondi2"], 
                    df_non_binding["ddg_gb_66_mbondi2"]) ,
                0.0 
                ] )

    ##### 

    tmp_df = pd.DataFrame.from_dict( 
                            { 
                                "Model" : ["gb5 mbondi2", "gb 66 mbondi2", "pb"] ,
                                "mean exp" : std_exp,
                                "RMSE exp" : mse_exp , 
                                 
                                "mean exp non binding" : std_exp_non_binding,
                                "RMSE exp non binding" : mse_exp_non_binding,
                                
                                "mean exp binding" : std_exp_binding,
                                "RMSE exp binding" : mse_exp_binding,
                                
                                "mean Pb" : std_pb ,
                                "RMSE Pb" : mse_pb ,
                                
                                "mean pb non binding" : std_pb_non_binding,
                                "RMSE pb non binding" : mse_pb_non_binding,
                                
                                "mean pb binding" : std_pb_binding,
                                "RMSE pb binding" : mse_pb_binding
                                
                            } )
    
        
    tmp_df.set_index('Model', inplace=True)
    return tmp_df




def radii_table_gen_mu_rmse(df) :
    exp_vals = df["ddg_exp"].to_numpy()
    mse_exp = np.array( [ pf.RMSE(exp_vals, df["ddg_gb_5_mbondi2"]) , 
                         pf.RMSE(exp_vals, df["ddg_gb_66_mbondi2"]),
                         pf.RMSE(exp_vals, df["ddg_gb_66_opt1"]),
                         pf.RMSE(exp_vals, df["ddg_gb_66_opt5d"]),
                         pf.RMSE(exp_vals, df["ddg_pb_mbondi2"]) ] )
    
    mu_exp = np.array( [ pf.mu(exp_vals, df["ddg_gb_5_mbondi2"]) , 
                        pf.mu(exp_vals, df["ddg_gb_66_mbondi2"]),
                        pf.mu(exp_vals, df["ddg_gb_66_opt1"]),
                        pf.mu(exp_vals, df["ddg_gb_66_opt5d"]),
                        pf.mu(exp_vals, df["ddg_pb_mbondi2"]) ] )
    
    
    mse_pb  = np.array( [ 
                         pf.RMSE(df["ddg_pb_mbondi2"], df["ddg_gb_5_mbondi2"]) , 
                         pf.RMSE(df["ddg_pb_mbondi2"], df["ddg_gb_66_mbondi2"]) , 
                         pf.RMSE(df["ddg_pb_mbondi2"], df["ddg_gb_66_opt1"]) , 
                         pf.RMSE(df["ddg_pb_mbondi2"], df["ddg_gb_66_opt5d"]) , 
                         0.0 ] )
    
    mu_pb  = np.array( [ 
                        pf.mu(df["ddg_pb_mbondi2"], df["ddg_gb_5_mbondi2"]) , 
                        pf.mu(df["ddg_pb_mbondi2"], df["ddg_gb_66_mbondi2"]) , 
                        pf.mu(df["ddg_pb_mbondi2"], df["ddg_gb_66_opt1"]) , 
                        pf.mu(df["ddg_pb_mbondi2"], df["ddg_gb_66_opt5d"]) , 
                        0.0 ] )

    
    
    df_binding = df[df["bd_site"] == True ].copy()
    mse_exp_binding = np.array( 
                    [ 
                    pf.RMSE(df_binding["ddg_exp"], 
                           df_binding["ddg_gb_5_mbondi2"]) , 
                    
                    pf.RMSE(df_binding["ddg_exp"], 
                           df_binding["ddg_gb_66_mbondi2"]),
                    pf.RMSE(df_binding["ddg_exp"], 
                           df_binding["ddg_gb_66_opt1"]),
                    pf.RMSE(df_binding["ddg_exp"], 
                           df_binding["ddg_gb_66_opt5d"]) ,
                    
                    pf.RMSE(df_binding["ddg_exp"], 
                           df_binding["ddg_pb_mbondi2"]) 
                    
                    ] )


    mu_exp_binding = np.array( 
                [ 
                pf.mu(df_binding["ddg_exp"], 
                        df_binding["ddg_gb_5_mbondi2"]) , 
                
                pf.mu(df_binding["ddg_exp"], 
                        df_binding["ddg_gb_66_mbondi2"]),
                pf.mu(df_binding["ddg_exp"], 
                        df_binding["ddg_gb_66_opt1"]),
                pf.mu(df_binding["ddg_exp"], 
                        df_binding["ddg_gb_66_opt5d"]) ,
                
                pf.mu(df_binding["ddg_exp"], 
                        df_binding["ddg_pb_mbondi2"]) 
                
                ] )

    mse_pb_binding  = np.array( 
                    [ 
                pf.RMSE(df_binding["ddg_pb_mbondi2"],
                       df_binding["ddg_gb_5_mbondi2"]) , 
                
                pf.RMSE(df_binding["ddg_pb_mbondi2"], 
                       df_binding["ddg_gb_66_mbondi2"]) ,
                
                pf.RMSE(df_binding["ddg_pb_mbondi2"], 
                           df_binding["ddg_gb_66_opt1"]),
                pf.RMSE(df_binding["ddg_pb_mbondi2"], 
                           df_binding["ddg_gb_66_opt5d"]) ,
                    0.0 
                    ] )
    
    mu_pb_binding  = np.array( 
                [ 
            pf.mu(df_binding["ddg_pb_mbondi2"],
                    df_binding["ddg_gb_5_mbondi2"]) , 
            
            pf.mu(df_binding["ddg_pb_mbondi2"], 
                    df_binding["ddg_gb_66_mbondi2"]) ,
            
            pf.mu(df_binding["ddg_pb_mbondi2"], 
                        df_binding["ddg_gb_66_opt1"]),
            pf.mu(df_binding["ddg_pb_mbondi2"], 
                        df_binding["ddg_gb_66_opt5d"]) ,
                0.0 
                ] )

    #non binding 
    df_non_binding = df[df["bd_site"] == False ].copy()
    mse_exp_non_binding = np.array( 
                    [ 
                    pf.RMSE(df_non_binding["ddg_exp"], 
                           df_non_binding["ddg_gb_5_mbondi2"]) , 
                    pf.RMSE(df_non_binding["ddg_exp"], 
                           df_non_binding["ddg_gb_66_mbondi2"]),
                    pf.RMSE(df_non_binding["ddg_exp"], 
                           df_non_binding["ddg_gb_66_opt1"]),
                    pf.RMSE(df_non_binding["ddg_exp"], 
                           df_non_binding["ddg_gb_66_opt5d"]) ,
                    
                    pf.RMSE(df_non_binding["ddg_exp"], 
                           df_non_binding["ddg_pb_mbondi2"]) 
                    
                    ] )

    mu_exp_non_binding = np.array( 
                [ 
                pf.mu(df_non_binding["ddg_exp"], 
                        df_non_binding["ddg_gb_5_mbondi2"]) , 
                pf.mu(df_non_binding["ddg_exp"], 
                        df_non_binding["ddg_gb_66_mbondi2"]),
                pf.mu(df_non_binding["ddg_exp"], 
                        df_non_binding["ddg_gb_66_opt1"]),
                pf.mu(df_non_binding["ddg_exp"], 
                        df_non_binding["ddg_gb_66_opt5d"]) ,
                
                pf.mu(df_non_binding["ddg_exp"], 
                        df_non_binding["ddg_pb_mbondi2"]) 
                
                ] )

    mse_pb_non_binding  = np.array( 
                    [ 
                pf.RMSE(df_non_binding["ddg_pb_mbondi2"],
                       df_non_binding["ddg_gb_5_mbondi2"]) ,
                 
                pf.RMSE(df_non_binding["ddg_pb_mbondi2"], 
                       df_non_binding["ddg_gb_66_mbondi2"]) ,
                
                pf.RMSE(df_non_binding["ddg_pb_mbondi2"], 
                        df_non_binding["ddg_gb_66_opt1"]),
                
                pf.RMSE(df_non_binding["ddg_pb_mbondi2"], 
                        df_non_binding["ddg_gb_66_opt5d"]) ,
                    
                    0.0 
                    ] )
    
    mu_pb_non_binding  = np.array( 
                [ 
            pf.mu(df_non_binding["ddg_pb_mbondi2"],
                    df_non_binding["ddg_gb_5_mbondi2"]) ,
                
            pf.mu(df_non_binding["ddg_pb_mbondi2"], 
                    df_non_binding["ddg_gb_66_mbondi2"]) ,
            
            pf.mu(df_non_binding["ddg_pb_mbondi2"], 
                    df_non_binding["ddg_gb_66_opt1"]),
            
            pf.mu(df_non_binding["ddg_pb_mbondi2"], 
                    df_non_binding["ddg_gb_66_opt5d"]) ,
                
                0.0 
                ] )
    
    ##### 

    tmp_df = pd.DataFrame.from_dict( 
                            { 
                                "Model" : ["gb5 mbondi2", 
                                           "gb 66 mbondi2",
                                           "gb 66 op1" ,
                                           "gb 66 op5d" ,
                                           "pb"] ,
                                "mean exp" : mu_exp ,
                                "RMSE exp" : mse_exp , 
                                
                                "mean exp non binding" : mu_exp_non_binding,
                                "RMSE exp non binding" : mse_exp_non_binding,
                                
                                "mean exp binding" : mu_exp_binding,
                                "RMSE exp binding" : mse_exp_binding,
                                
                                "mean Pb" : mu_pb ,
                                "RMSE Pb" : mse_pb ,
                                
                                "mean pb non binding" : mu_pb_non_binding,
                                "RMSE pb non binding" : mse_pb_non_binding,
                                
                                "mean pb binding" : mu_pb_binding,
                                "RMSE pb binding" : mse_pb_binding
                            } )
    
        
    tmp_df.set_index('Model', inplace=True)
    return tmp_df



def row_helper_exp(df, func ) :
       func_out = np.array( [ func(df["ddg_exp"], df["ddg_gb_5_mbondi2"]) , 
                            func(df["ddg_exp"], df["ddg_gb_66_mbondi2"]),
                            func(df["ddg_exp"], df["ddg_pb_mbondi2"]) ] )
       return func_out

def row_helper_pb(df,func ) :
       func_out  = np.array( [ 
                     func(df["ddg_pb_mbondi2"], df["ddg_gb_5_mbondi2"]) , 
                     func(df["ddg_pb_mbondi2"], df["ddg_gb_66_mbondi2"]) , 
                            0.0 ] )
       return func_out

def corr_helper_exp(df):
	corr = df[["ddg_exp" ,"ddg_gb_5_mbondi2","ddg_gb_66_mbondi2", "ddg_pb_mbondi2" ]].corr()
	return corr["ddg_exp"].to_numpy()[1:]
def corr_helper_pb(df):
	corr = df[["ddg_exp" ,"ddg_gb_5_mbondi2","ddg_gb_66_mbondi2", "ddg_pb_mbondi2" ]].corr()
	return corr["ddg_pb_mbondi2"].to_numpy()[1:]
def wrost_20_rmse(arr1, arr2) :
	error_array = (arr1 - arr2)**2
	threshold = np.percentile(error_array, 80)
	error_worst_20 = error_array[error_array >= threshold]
	return np.sqrt(np.mean(error_worst_20))

def table_gen_rmse_corr_bounds(df, split_by_pb=False) :
	#exp_vals = df["ddg_exp"].to_numpy()
	mse_exp = row_helper_exp(df, pf.RMSE)
	mse_pb  = row_helper_pb(df, pf.RMSE)
	wrost_20_rmse_exp = row_helper_exp(df, wrost_20_rmse)
	full_corr_exp= corr_helper_exp(df)
 
	mean_exp = row_helper_exp(df, pf.mu)
	mean_pb = row_helper_pb(df, pf.mu)
	full_corr_pb = corr_helper_pb(df)
	wrost_20_rmse_pb = row_helper_pb(df, wrost_20_rmse)
	

 
 
	df_binding = df[df["bd_site"] == True ].copy()
	mse_exp_binding = row_helper_exp(df_binding, pf.RMSE)
	mse_pb_binding  = row_helper_pb(df_binding, pf.RMSE)
	mean_exp_binding = row_helper_exp(df_binding, pf.mu)
	mean_pb_binding = row_helper_pb(df_binding, pf.mu)
	wrost_20_rmse_exp_binding = row_helper_exp(df_binding, wrost_20_rmse)
	wrost_20_rmse_pb_binding = row_helper_pb(df_binding, wrost_20_rmse)
	binding_corr_exp = corr_helper_exp(df_binding)
	binding_corr_pb = corr_helper_pb(df_binding)
	#non binding 
	df_non_binding = df[df["bd_site"] == False ].copy() 
	mse_exp_non_binding = row_helper_exp(df_non_binding, pf.RMSE)
	mse_pb_non_binding  = row_helper_pb(df_non_binding, pf.RMSE)
	mean_exp_non_binding = row_helper_exp(df_non_binding, pf.mu)
	mean_pb_non_binding = row_helper_pb(df_non_binding, pf.mu)
	wrost_20_rmse_exp_non_binding = row_helper_exp(df_non_binding, wrost_20_rmse)
	wrost_20_rmse_pb_non_binding = row_helper_pb(df_non_binding, wrost_20_rmse)
	non_binding_corr_exp = corr_helper_exp(df_non_binding)
	non_binding_corr_pb = corr_helper_pb(df_non_binding)

	##### 
       
	if split_by_pb ==False:
		output_df = pd.DataFrame.from_dict({ "Model" : ["GB5 Mbondi2", "GB66 Mbondi2", "PB"] ,
				r"RMSE $\Delta \Delta  \Delta G_{\text{bind,exp}}$ BD " : mse_exp ,
				r"mean $\Delta \Delta  \Delta G_{\text{bind,exp}}$ BD" : mean_exp,
				r"corr $\Delta \Delta  \Delta G_{\text{bind,exp}}$ BD" : full_corr_exp,
				r"RMSE of worst 20\% $\Delta \Delta  \Delta G_{\text{bind,exp}}$ BD": wrost_20_rmse_exp,
				
				r"RMSE $\Delta \Delta  \Delta G_{\text{bind,PB}}$ BD" : mse_pb ,
				r"mean $\Delta \Delta  \Delta G_{\text{bind,PB}}$ BD" : mean_pb ,
				r"corr $\Delta \Delta  \Delta G_{\text{bind,PB}}$ BD" : full_corr_pb,
				r"RMSE of worst 20\% $\Delta \Delta  \Delta G_{\text{bind,PB}}$ BD": wrost_20_rmse_pb,
				
				r"RMSE  $\Delta \Delta  \Delta G_{\text{bind,exp}}$ BS" : mse_exp_binding,
				r"mean  $\Delta \Delta  \Delta G_{\text{bind,exp}}$ BS" : mean_exp_binding,
				r"corr  $\Delta \Delta  \Delta G_{\text{bind,exp}}$ BS" : binding_corr_exp,
				r"RMSE of worst 20\%  $\Delta \Delta  \Delta G_{\text{bind,exp}}$ BS": wrost_20_rmse_exp_binding,
				
				r"RMSE $\Delta \Delta  \Delta G_{\text{bind,PB}}$ BS" : mse_pb_binding,
				r"mean $\Delta \Delta  \Delta G_{\text{bind,PB}}$ BS" : mean_pb_binding,
				r"corr $\Delta \Delta  \Delta G_{\text{bind,PB}}$ BS" : binding_corr_pb,
				r"RMSE of worst 20\% $\Delta \Delta  \Delta G_{\text{bind,PB}}$ BS": wrost_20_rmse_pb_binding,
				
				r"RMSE $\Delta \Delta  \Delta G_{\text{bind,exp}}$ N-BS" : mse_exp_non_binding,
				r"mean $\Delta \Delta  \Delta G_{\text{bind,exp}}$ N-BS" : mean_exp_non_binding,
				r"corr $\Delta \Delta  \Delta G_{\text{bind,exp}}$ N-BS" : non_binding_corr_exp,
				r"RMSE of worst 20\% $\Delta \Delta  \Delta G_{\text{bind,exp}}$ N-BS": wrost_20_rmse_exp_non_binding,
	
				r"RMSE $\Delta \Delta  \Delta G_{\text{bind,PB}}$ N-BS" : mse_pb_non_binding,
				r"mean $\Delta \Delta  \Delta G_{\text{bind,PB}}$  N-BS" : mean_pb_non_binding,
				r"corr $\Delta \Delta  \Delta G_{\text{bind,PB}}$  N-BS" : non_binding_corr_pb,
				r"RMSE of worst 20\% $\Delta \Delta  \Delta G_{\text{bind,PB}}$ N-BS": wrost_20_rmse_pb_non_binding
		
				})
	
		output_df.set_index('Model', inplace=True)
		return output_df
	else:
		output_df_exp = pd.DataFrame.from_dict({ "Model" : ["GB5 Mbondi2", "GB66 Mbondi2", "PB"],
				r"RMSE $\Delta \Delta  \Delta G_{\text{bind,exp}}$ BD " : mse_exp ,
				r"mean $\Delta \Delta  \Delta G_{\text{bind,exp}}$ BD" : mean_exp,
				r"corr $\Delta \Delta  \Delta G_{\text{bind,exp}}$ BD" : full_corr_exp,
				r"RMSE of worst 20\% $\Delta \Delta  \Delta G_{\text{bind,exp}}$ BD": wrost_20_rmse_exp,


				r"RMSE $\Delta \Delta  \Delta G_{\text{bind,exp}}$ N-BS" : mse_exp_non_binding,
				r"mean $\Delta \Delta  \Delta G_{\text{bind,exp}}$ N-BS" : mean_exp_non_binding,
				r"corr $\Delta \Delta  \Delta G_{\text{bind,exp}}$ N-BS" : non_binding_corr_exp,
				r"RMSE of worst 20\% $\Delta \Delta  \Delta G_{\text{bind,exp}}$ N-BS": wrost_20_rmse_exp_non_binding,
    
                            				
			       r"RMSE  $\Delta \Delta  \Delta G_{\text{bind,exp}}$ BS" : mse_exp_binding,
				r"mean  $\Delta \Delta  \Delta G_{\text{bind,exp}}$ BS" : mean_exp_binding,
				r"corr  $\Delta \Delta  \Delta G_{\text{bind,exp}}$ BS" : binding_corr_exp,
				r"RMSE of worst 20\%  $\Delta \Delta  \Delta G_{\text{bind,exp}}$ BS": wrost_20_rmse_exp_binding
				})
		output_df_exp.set_index('Model', inplace=True)

		output_df_pb = pd.DataFrame.from_dict({ "Model" : ["GB5 Mbondi2", "GB66 Mbondi2", "PB"] ,
				r"RMSE $\Delta \Delta  \Delta G_{\text{bind,PB}}$ BD" : mse_pb ,
				r"mean $\Delta \Delta  \Delta G_{\text{bind,PB}}$ BD" : mean_pb ,
				r"corr $\Delta \Delta  \Delta G_{\text{bind,PB}}$ BD" : full_corr_pb,
				r"RMSE of worst 20\% $\Delta \Delta  \Delta G_{\text{bind,PB}}$ BD": wrost_20_rmse_pb,
				
                            r"RMSE $\Delta \Delta  \Delta G_{\text{bind,PB}}$ N-BS" : mse_pb_non_binding,
				r"mean $\Delta \Delta  \Delta G_{\text{bind,PB}}$  N-BS" : mean_pb_non_binding,
				r"corr $\Delta \Delta  \Delta G_{\text{bind,PB}}$  N-BS" : non_binding_corr_pb,
				r"RMSE of worst 20\% $\Delta \Delta  \Delta G_{\text{bind,PB}}$ N-BS": wrost_20_rmse_pb_non_binding,
				
				r"RMSE $\Delta \Delta  \Delta G_{\text{bind,PB}}$ BS" : mse_pb_binding,
				r"mean $\Delta \Delta  \Delta G_{\text{bind,PB}}$ BS" : mean_pb_binding,
				r"corr $\Delta \Delta  \Delta G_{\text{bind,PB}}$ BS" : binding_corr_pb,
				r"RMSE of worst 20\% $\Delta \Delta  \Delta G_{\text{bind,PB}}$ BS": wrost_20_rmse_pb_binding
				

				})
		output_df_pb.set_index('Model', inplace=True)
		return output_df_exp, output_df_pb
	

		
	

def row_helper_radii_exp(df, func ) :
    func_out = np.array( 
                    [ 
                    func(df["ddg_exp"], 
                           df["ddg_gb_5_mbondi2"]) , 
                    
                    func(df["ddg_exp"], 
                           df["ddg_gb_66_mbondi2"]),
                    func(df["ddg_exp"], 
                           df["ddg_gb_66_opt1"]),
                    func(df["ddg_exp"], 
                           df["ddg_gb_66_opt5d"]) ,
                    
                    func(df["ddg_exp"], 
                           df["ddg_pb_mbondi2"]) 
                    
                    ] )
    return func_out


def row_helper_radii_pb(df,func ) :
	func_out  = np.array( 
					[ 
				func(df["ddg_pb_mbondi2"],
						df["ddg_gb_5_mbondi2"]) , 
				
				func(df["ddg_pb_mbondi2"], 
						df["ddg_gb_66_mbondi2"]) ,
				
				func(df["ddg_pb_mbondi2"], 
							df["ddg_gb_66_opt1"]),
				func(df["ddg_pb_mbondi2"], 
							df["ddg_gb_66_opt5d"]) ,
					0.0 
					] )
	return func_out


def corr_func(*args) :
       return np.corrcoef(*args)[0,1] #get off diag
     
def full_radii_table_gen_mu_rmse(df, split_by_pb=False) :
       corc_coef_func = corr_func
       mse_exp = row_helper_radii_exp(df, pf.RMSE)
       mu_exp = row_helper_radii_exp(df, pf.mu)
       rmse_worst_20_exp = row_helper_radii_exp(df, wrost_20_rmse)
       corr_exp = row_helper_radii_exp(df,corc_coef_func)



       mse_pb  = row_helper_radii_pb(df, pf.RMSE)
       mu_pb  = row_helper_radii_pb(df, pf.mu)
       rmse_worst_20_pb = row_helper_radii_pb(df, wrost_20_rmse)
       corr_pb = row_helper_radii_pb(df, corc_coef_func)


       df_binding = df[df["bd_site"] == True ].copy()
       
       mse_exp_binding = row_helper_radii_exp(df_binding, pf.RMSE)
       mu_exp_binding = row_helper_radii_exp(df_binding, pf.mu)
       rmse_worst_20_exp_binding = row_helper_radii_exp(df_binding, wrost_20_rmse)
       corr_exp_binding = row_helper_radii_exp(df_binding,corc_coef_func)
       
       
       
       mse_pb_binding  = row_helper_radii_pb(df_binding, pf.RMSE)
       mu_pb_binding  = row_helper_radii_pb(df_binding, pf.mu)
       rmse_worst_20_pb_binding = row_helper_radii_pb(df_binding, wrost_20_rmse)
       corr_pb_binding = row_helper_radii_pb(df_binding,corc_coef_func)
       
       
       #non binding 
       df_non_binding = df[df["bd_site"] == False ].copy()
       
       
       mse_exp_non_binding = row_helper_radii_exp(df_non_binding, pf.RMSE)
       mu_exp_non_binding = row_helper_radii_exp(df_non_binding, pf.mu)
       rmse_worst_20_exp_non_binding = row_helper_radii_exp(df_non_binding, wrost_20_rmse)
       corr_exp_non_binding = row_helper_radii_exp(df_non_binding, corc_coef_func)


       mse_pb_non_binding  = row_helper_radii_pb(df_non_binding, pf.RMSE)
       mu_pb_non_binding  = row_helper_radii_pb(df_non_binding, pf.mu)
       rmse_worst_20_pb_non_binding = row_helper_radii_pb(df_non_binding, wrost_20_rmse)
       corr_pb_non_binding = row_helper_radii_pb(df_non_binding, corc_coef_func)

       ##### 
       
       pb_str = r"$\Delta \Delta  \Delta G_{\text{bind,PB}}$"
       exp_str = r"$\Delta \Delta  \Delta G_{\text{bind,exp}}$"
       if split_by_pb == False:
              tmp_df = pd.DataFrame.from_dict({
                     "Model" :["GB5 Mbondi2",
                                          "GB66 Mbondi2",
                                          "GB66  OPT1" ,
                                          "GB66  OPT5D" ,
                                          "PB"]  ,
                     
                     f"RMSE {exp_str} BD" : mse_exp ,
                     f"mean {exp_str} BD" : mu_exp ,
                     f"PCC {exp_str} BD" : corr_exp, 
                     f"RMSE of worst 20\% {exp_str} BD": rmse_worst_20_exp,

                     f"RMSE {pb_str} BD" : mse_pb ,
                     f"mean {pb_str} BD" : mu_pb ,
                     f"PCC {pb_str} BD" : corr_pb, 
                     f"RMSE of worst 20\% {pb_str} BD": rmse_worst_20_pb,

                     
                     f"RMSE {exp_str} N-BS" : mse_exp_non_binding,
                     f"mean {exp_str} N-BS" : mu_exp_non_binding,
                     f"PCC {exp_str} N-BS" : corr_exp_non_binding, 
                     f"RMSE of worst 20\% {exp_str} N-BS": rmse_worst_20_exp_non_binding,
                     
                     f"RMSE {pb_str} N-BS" : mse_pb_non_binding,
                     f"mean {pb_str} N-BS" : mu_pb_non_binding,
                     f"PCC {pb_str} N-BS" : corr_pb_non_binding, 
                     f"RMSE of worst 20\% {pb_str} N-BS": rmse_worst_20_pb_non_binding,
                     
                     f"RMSE {exp_str} BS" : mse_exp_binding,
                     f"mean {exp_str} BS" : mu_exp_binding,
                     f"PCC {exp_str} BS" : corr_exp_binding, 
                     f"RMSE of worst 20\% {exp_str} BS": rmse_worst_20_exp_binding,

                     f"mean {pb_str} BS" : mu_pb_binding,
                     f"RMSE {pb_str} BS" : mse_pb_binding,
                     f"PCC {pb_str} BS" : corr_pb_binding, 
                     f"RMSE of worst 20\% {pb_str} BS": rmse_worst_20_pb_binding

                     

              })
              tmp_df.set_index('Model', inplace=True)
              return tmp_df
       else:
              exp_ref_df = pd.DataFrame.from_dict({
                     "Model" :["GB5 Mbondi2",
                                          "GB66 Mbondi2",
                                          "GB66  OPT1" ,
                                          "GB66  OPT5D" ,
                                          "PB"]  ,
                     f"RMSE {exp_str} BD" : mse_exp ,
                     f"mean {exp_str} BD" : mu_exp ,
                     f"PCC {exp_str} BD" : corr_exp, 
                     f"RMSE of worst 20\% {exp_str} BD": rmse_worst_20_exp,


                     f"RMSE {exp_str} N-BS" : mse_exp_non_binding,
                     f"mean {exp_str} N-BS" : mu_exp_non_binding,
                     f"PCC {exp_str} N-BS" : corr_exp_non_binding, 
                     f"RMSE of worst 20\% {exp_str} N-BS": rmse_worst_20_exp_non_binding,

                     f"RMSE {exp_str} BS" : mse_exp_binding,
                     f"mean {exp_str} BS" : mu_exp_binding,
                     f"PCC {exp_str} BS" : corr_exp_binding, 
                     f"RMSE of worst 20\% {exp_str} BS": rmse_worst_20_exp_binding,
              })
              exp_ref_df.set_index('Model', inplace=True)
              pb_ref_df = pd.DataFrame.from_dict({
                     "Model" : ["GB5 Mbondi2",
                                          "GB66 Mbondi2",
                                          "GB66  OPT1" ,
                                          "GB66  OPT5D" ,
                                          "PB"] ,
                     f"RMSE {pb_str} BD" : mse_pb ,
                     f"mean {pb_str} BD" : mu_pb ,
                     f"PCC {pb_str} BD" : corr_pb, 
                     f"RMSE of worst 20\% {pb_str} BD": rmse_worst_20_pb,

                     f"RMSE {pb_str} N-BS" : mse_pb_non_binding,
                     f"mean {pb_str} N-BS" : mu_pb_non_binding,
                     f"PCC {pb_str} N-BS" : corr_pb_non_binding, 
                     f"RMSE of worst 20\% {pb_str} N-BS": rmse_worst_20_pb_non_binding,
                     
                     f"mean {pb_str} BS" : mu_pb_binding,
                     f"RMSE {pb_str} BS" : mse_pb_binding,
                     f"PCC {pb_str} BS" : corr_pb_binding, 
                     f"RMSE of worst 20\% {pb_str} BS": rmse_worst_20_pb_binding

              })
              pb_ref_df.set_index('Model', inplace=True)
              return exp_ref_df, pb_ref_df





def calculate_features(df: pd.DataFrame, 
                       columns: list,
                       ):
	"""
	Calculates various statistical features for specified columns in a DataFrame.

	Args:
		df (pd.DataFrame): The input DataFrame.
		columns (list): A list of column names to analyze.
              

	Returns:
		pd.DataFrame: A DataFrame containing the calculated statistics
						for each specified column.
	"""
	results = {}
	for col in columns:
		# positive_data = df[df[col] > 0][col]
		# relative_mass_positive = len(positive_data) / len(df) if len(df) > 0 else 0
		positive_data = df[df[col] > 0][col]
		pos_mass = np.sum(positive_data)
		total_mass = np.sum(np.abs(df[col]))
		relative_mass_positive = pos_mass / total_mass 
  
  
		mean_all = df[col].mean()
		median_all = df[col].median()
		mode_all = df[col].mode().iloc[0] if not df[col].mode().empty else None
		pskew_all = skew(df[col])
		std_all = df[col].std()

		results[col] = {
			"PPRM": relative_mass_positive,
			"mean": mean_all,
			"median": median_all,
			"mode": mode_all,
			"p_skewness": pskew_all,
			"std": std_all,
		}
		
		
	df_truncated = df[df["bd_site"] == False]
	for col in columns:
		positive_data = df_truncated[df_truncated[col] > 0][col]
		pos_mass = np.sum(positive_data)
		total_mass = np.sum(np.abs(df_truncated[col]))
		relative_mass_positive = pos_mass / total_mass 
  
		# relative_mass_positive = len(positive_data) / len(df_truncated) if len(df_truncated) > 0 else 0

		mean_all = df_truncated[col].mean()
		median_all = df_truncated[col].median()
		mode_all = df_truncated[col].mode().iloc[0] if not df_truncated[col].mode().empty else None
		pskew_all = skew(df_truncated[col])
		std_all = df_truncated[col].std()

		results[f"{col}_NBS"] = {
			"PPRM": relative_mass_positive,
			"mean": mean_all,
			"median": median_all,
			"mode": mode_all,
			"p_skewness": pskew_all,
			"std": std_all,
		}

	df_truncated = df[df["bd_site"] == True].copy()
	for col in columns:
		# positive_data = df_truncated[df_truncated[col] > 0][col]
		# relative_mass_positive = len(positive_data) / len(df_truncated) if len(df_truncated) > 0 else 0

		positive_data = df_truncated[df_truncated[col] > 0][col]
		pos_mass = np.sum(positive_data)
		total_mass = np.sum(np.abs(df_truncated[col]))
		relative_mass_positive = pos_mass / total_mass 

		mean_all = df_truncated[col].mean()
		median_all = df_truncated[col].median()
		mode_all = df_truncated[col].mode().iloc[0] if not df_truncated[col].mode().empty else None
		pskew_all = skew(df_truncated[col])
		std_all = df_truncated[col].std()

		results[f"{col}_BS"] = {
			"PPRM": relative_mass_positive,
			"mean": mean_all,
			"median": median_all,
			"mode": mode_all,
			"p_skewness": pskew_all,
			"std": std_all,
		}

	return pd.DataFrame.from_dict(results, orient='index')





def _df_to_latex_cust(df, bold=True):
       if bold: 
              latex_output = (
                     df.style
                     .format("{:.2f}") # Apply float formatting using Styler.format()
                     .highlight_min(axis=1, props='textbf:--rwrap;') # Or axis=0 for min in columns
                     .to_latex(
                            hrules=True # Optional: Adds \toprule, \midrule, \bottomrule
                     )
              )
              print(latex_output)
       else:
              latex_output = (
                     df.style
                     .format("{:.2f}") # Apply float formatting using Styler.format()
                     .to_latex(
                            hrules=True # Optional: Adds \toprule, \midrule, \bottomrule
                     )
              )
              print(latex_output)
		

		
       
        
def neut_table_gen(df) :
    exp_vals = df["ddg_exp"].to_numpy()
    mse_exp = np.array( [ pf.RMSE(exp_vals, df["ddg_gb_5_mbondi2"]) , 
                         pf.RMSE(exp_vals, df["ddg_gb_66_mbondi2"]),
                         
                         
                         pf.RMSE(exp_vals, df["ddg_gb_5_Nmbondi2"]) , 
                         pf.RMSE(exp_vals, df["ddg_gb_66_Nmbondi2"]),
                         
                         ] )


    df_binding = df[df["bd_site"] == True ].copy()
    #bad code, but w/e 
    mse_exp_binding = np.array( 
                    [ 
                    pf.RMSE(df_binding["ddg_exp"], 
                           df_binding["ddg_gb_5_mbondi2"]) , 
                    pf.RMSE(df_binding["ddg_exp"], 
                           df_binding["ddg_gb_66_mbondi2"]),
                    
                    
                    pf.RMSE(df_binding["ddg_exp"], 
                           df_binding["ddg_gb_5_Nmbondi2"]) , 
                    pf.RMSE(df_binding["ddg_exp"], 
                           df_binding["ddg_gb_66_Nmbondi2"]),
                    
                    ] )



    #non binding 
    df_non_binding = df[df["bd_site"] == False ].copy()
    mse_exp_non_binding = np.array( 
                    [ 
                    pf.RMSE(df_non_binding["ddg_exp"], 
                           df_non_binding["ddg_gb_5_mbondi2"]) , 
                    pf.RMSE(df_non_binding["ddg_exp"], 
                           df_non_binding["ddg_gb_66_mbondi2"]),
                    
                    pf.RMSE(df_non_binding["ddg_exp"], 
                           df_non_binding["ddg_gb_5_Nmbondi2"]) , 
                    pf.RMSE(df_non_binding["ddg_exp"], 
                           df_non_binding["ddg_gb_66_Nmbondi2"]),
                    ] )
    
    
    ###   
    tmp_df = pd.DataFrame.from_dict( 
                            { 
                                "Model" : ["gb5 mbondi2", "gb 66 mbondi2", "gb5 Nmbondi2", "gb 66 Nmbondi2",] ,
                                "RMSE exp" : mse_exp , 
                                "RMSE exp N-BS" : mse_exp_non_binding,
                                "RMSE exp binding" : mse_exp_binding,
                            } )
    
        
    tmp_df.set_index('Model', inplace=True)
    return tmp_df