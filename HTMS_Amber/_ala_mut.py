import os
import re

import _defaults
import _utils


def mutations(pdb_data, name_from, name_to, idx) -> list:
    '''
    pdb_data    : pdb file we want to mutate
    name_from   : three letter name for the initial amino acid 
    name_to     : three letter name for the final amino acid
    idx         : index of the mutation 
    
    returns     : mutated pdb as a list 
    '''
    counter = 0 
    records = ('ATOM', 'HETATM', 'TER', 'ANISOU')
    mutated_ligand = [] 
    
    for line in pdb_data:
        if line.startswith(records):
            
            
            if (line[17:20].strip() in name_from) and (line[22:26].strip() == idx):
                #add cases here for longer ones... 
                #add case for GLN which has AGLN and BGLN
                if (counter <= 4) and (line[16] == "A" or line[16]== " ") : #count for ALA:  
                    new_line= line[:16] + name_to.rjust(4) + line[20:]
                    mutated_ligand.append(new_line)
                    counter = counter + 1 
                    continue
                else :
                    continue
        mutated_ligand.append(line)

    return mutated_ligand


def pdb_split(pdb_data, option) -> list:
    '''
    pdb_data    : pdb file we want to split
    option      : option will determine if we want to get the receptor or the ligand,
    receptor =0, ligand =1 
     
    returns     : split as a list 
    '''
    #ignore HET and other line starts: 
    #NOTE: no need to use the no_HET source files, this will strip the files of the HET
    ter_state = 0 
    records = ('ATOM', 'ANISOU', 'TER')
    data = []
    
    for line in pdb_data:
        if line.startswith(records):
            if (option == 0) and (ter_state==0)  : 
                data.append(line)
                
                if line.startswith('TER') :
                    return data #break once we get to first Ter as option 0
            
            #need to check for Ter after store line starting with Ter due to structure 
            #of pdb files, TER line belongs to structure. 
            if line.startswith('TER') and (ter_state==0) :
                ter_state =1
                continue 
            if (option == 1 ) and (ter_state==1): 
                data.append(line)
                
    return data


def split_and_mut(pdbfh, pdbfh_base_name, name_from, name_to, idx, naming_conv) :
    with open(pdbfh, "r") as f :
        pdb_data = f.readlines()
    #splits
    struct_pdb_data = pdb_split(pdb_data, 0 )
    file_handle_structure = pdbfh_base_name + "_recpt.pdb"
    with open(file_handle_structure, "w+") as pdb_file : 
        for line in struct_pdb_data : 
            pdb_file.write(f"{line}")
        pdb_file.close()
    ligand_pdb = pdb_split(pdb_data, 1 )
    file_handle_ligand = pdbfh_base_name + "_ligand.pdb"
    with open(file_handle_ligand, "w+") as pdb_file : 
        for line in ligand_pdb : 
            pdb_file.write(f"{line}")
        pdb_file.close()
    #mutations:
    #ligand_mutation
    mutation_pdb_data = mutations(ligand_pdb, name_from, name_to, idx)
    file_handle_mut_base = pdbfh_base_name +"_" + naming_conv
    file_handle_mut = file_handle_mut_base + "_ligand.pdb"
    with open(file_handle_mut, "w+") as pdb_file : 
        for line in mutation_pdb_data : 
            pdb_file.write(f"{line}")
        pdb_file.close()
    #full file mutation
    mutation_pdb_data_all = mutations(pdb_data, name_from, name_to, idx) 
    file_handle_mut_all = file_handle_mut_base + ".pdb"
    with open(file_handle_mut_all, "w+") as pdb_file : 
        for line in mutation_pdb_data_all : 
            pdb_file.write(f"{line}")
        pdb_file.close()
    return  file_handle_mut_base

def _extract_mut_info(s):
    pattern = r'([A-Za-z])(\d+)([A-Za-z])'
    match = re.search(pattern, s)
    if match:
        name_from_char = str(match.group(1))
        idx = str(match.group(2))
        name_to_char = str(match.group(3))
    else:
        raise ValueError(f"Could not extract residue number from string: {s}")
    return name_from_char, idx, name_to_char
 
def general_method(input_dict : dict, pdbfh, 
                   pdbfh_base_name, mutation, just_build, amber_source) : 
    """
    general process: 
    input_dict: from input_args_check
    mut_num : range(len(input_dict["MUTATIONS"]))
    """
    
    name_from_char, idx, name_to_char = _extract_mut_info(mutation)
    #dict for char to code conversion 
    amino_acid_dict =_utils.amino_acids 
    #convert to 3 letter code
    name_from   = amino_acid_dict[name_from_char]
    name_to     = amino_acid_dict[name_to_char]
    
    #get cwd
    cwd = os.getcwd()
    #get the three letter code as a str e.g. E484A
    naming_conv = name_from_char+ idx + name_to_char 
    #make a directory named: base_name_naming-conv_dir
    #string for dir name
    dir_name = pdbfh_base_name + "_" + naming_conv + "_dir"
    #path to dir
    dir_name_path =os.path.join(".", dir_name)
    #path for other files
  
    os.system(f"mkdir {dir_name_path}")
    
    dir_name_path_full = dir_name_path + "/"
    #new name for pdb in the dir
    pdbfh_in_dir = dir_name_path_full + pdbfh
    
    #copy the base pdb into new dir
    os.system(f"cp {pdbfh} {pdbfh_in_dir}")
    #update base name to the file in subdir
    #pdbfh_base_name = pdbfh_base_name_in_dir
    os.chdir(dir_name_path) #note the change back use chdir("..")
    
    ##TODO make into smaller functions 
    file_handle_mut_base = split_and_mut(pdbfh,
                                         pdbfh_base_name, 
                                         name_from, 
                                         name_to, 
                                         idx,
                                         naming_conv)
    
    ######################### tleap gen ######################################
    
    
    tleap_file_name = _defaults.tleap_in_gen( input_dict,
                                    pdbfh_base_name, 
                                    file_handle_mut_base)
    

    ########################## MMPBSA sh file gen #################
    run_MMPBSA_sh_name = _defaults.mmbpsa_sh_gen(input_dict , 
                                       pdbfh_base_name , 
                                       file_handle_mut_base,
                                       cwd, amber_source)

    if just_build:
        pass
    else:
        os.system(f"sbatch {run_MMPBSA_sh_name}")
    os.chdir("..")
    return 

