import os


def tleap_gen(pdbfh_base_name,file_handle_mut_all ) -> list:
    '''
    pdbfh_base_name     : base name of the pdb file
    file_handle_mut_all : base bame of the mutated file 
    returns             : tleap file as a list 
    '''

    tleap_wild_in = [f"source leaprc.protein.ff19SB",
        f"source leaprc.water.opc",
        f"set default PBRadii mbondi2\n",
        f"com = loadpdb {pdbfh_base_name}.pdb"  ,
        f"ligand = loadpdb {pdbfh_base_name}_ligand.pdb" ,
        f"rcp = loadpdb {pdbfh_base_name}_recpt.pdb\n",
        f"saveamberparm com {pdbfh_base_name}.prmtop {pdbfh_base_name}.inpcrd",
        f"saveamberparm ligand {pdbfh_base_name}_ligand.prmtop {pdbfh_base_name}_ligand.inpcrd",
        f"saveamberparm rcp {pdbfh_base_name}_recpt.prmtop {pdbfh_base_name}_recpt.inpcrd",
        f"com_md = loadpdb {pdbfh_base_name}.pdb"  ,
        f"solvatebox com_md OPCBOX 12.0",
        f"addIons2 com_md Na+ 0",
        f"saveamberparm com_md {pdbfh_base_name}_solvated.prmtop {pdbfh_base_name}_solvated.inpcrd\n",
        f"com_mut = loadpdb {file_handle_mut_all}.pdb",
        f"ligand_mut = loadpdb {file_handle_mut_all}_ligand.pdb\n",
        f"saveamberparm com_mut {file_handle_mut_all}.prmtop {file_handle_mut_all}.inpcrd",
        f"saveamberparm ligand_mut {file_handle_mut_all}_ligand.prmtop {file_handle_mut_all}_ligand.inpcrd",
        f"quit"]
    return tleap_wild_in

def tleap_in_gen( input_dict,  pdbfh_base_name, file_handle_mut_base ): 
    if input_dict["LEAP.IN_PATH" ] == [] :
        tleap_mut_in = tleap_gen(pdbfh_base_name, file_handle_mut_base)
        with open("tleap_mut.in", "w+") as tleap : 
            for line in tleap_mut_in : 
                tleap.write(f"{line}\n")
            tleap.close()
        os.system(f"dos2unix tleap_mut.in") #not sure if needed. 
        tleap_file_name ="tleap_mut.in"
    else : 
        tleap_file_name =input_dict["LEAP.IN_PATH" ]
    return  tleap_file_name

def mut_bash( pdbfh_base_name, 
             file_handle_mut_all, 
             cwd, 
             amber_source,
             tleap_file_name : str ='tleap_mut.in') -> list:
    '''
    pdbfh_base_name     : base name of the pdb file
    file_handle_mut_all : base bame of the mutated file 
    cwd                 : cwd where script was called, this is parent dir afer chdir call
    returns             : .sh file as a list 
    '''
    
    mut_bash_sh = [f"#!/bin/bash",
        f"#SBATCH --job-name=run_{pdbfh_base_name}_mut",
        f"#SBATCH --partition=cpu",
        f"#SBATCH --ntasks=4",
        f"#SBATCH --cpus-per-task=1",
        f"#SBATCH --mem=10000",
        f"#SBATCH --output=run_mmpbsa_{pdbfh_base_name}.out",
        f"#SBATCH --error=run_mmpbsa_{pdbfh_base_name}.error",
        f"#SBATCH --time=72:00:00",
        f'''echo "Loading modules..."'''  ,  
        f"module load amber " ,
        f"source {amber_source}",
        f"",
        f"tleap -s -f {tleap_file_name} > tleap_mut.out"
        f"",
        f"""$AMBERHOME/bin/MMPBSA.py -O -i \
{cwd}/mmpbsa.in -o \
FINAL_RESULTS_MMPBSA_tleap_{file_handle_mut_all}.dat\
 -sp {pdbfh_base_name}_solvated.prmtop\
 -cp {pdbfh_base_name}.prmtop\
 -rp {pdbfh_base_name}_recpt.prmtop\
 -lp {pdbfh_base_name}_ligand.prmtop\
 -y {cwd}/*.mdcrd\
 -mc {file_handle_mut_all}.prmtop\
 -ml {file_handle_mut_all}_ligand.prmtop"""
        ]
    return mut_bash_sh

    # f"{cwd}/change_radii_to_opt.py {pdbfh_base_name}_solvated.prmtop",
    #     f"{cwd}/change_radii_to_opt.py {pdbfh_base_name}.prmtop",
    #     f"{cwd}/change_radii_to_opt.py {pdbfh_base_name}_ligand.prmtop",
    #     f"{cwd}/change_radii_to_opt.py {file_handle_mut_all}.prmtop", 
    #     f"{cwd}/change_radii_to_opt.py {file_handle_mut_all}_ligand.prmtop",
    #     f"{cwd}/change_radii_to_opt.py {pdbfh_base_name}_recpt.prmtop",

def mmbpsa_sh_gen(input_dict , pdbfh_base_name,
                  file_handle_mut_base, cwd,
                  amber_source : str = "amber22/amber.sh"):
    if input_dict["MMPBSA.SH_PATH"] == []:
        mut_bash_file = mut_bash(pdbfh_base_name, 
                                 file_handle_mut_base, 
                                 cwd, amber_source)
        
        with open("run_MMPBSA.sh", "w+") as mut_bash_sh : 
            for line in mut_bash_file : 
                mut_bash_sh.write(f"{line}\n")
            mut_bash_sh.close()
            run_MMPBSA_sh_name = "run_MMPBSA.sh"
    else :
        run_MMPBSA_sh_name = input_dict["MMPBSA.SH_PATH"]
    return run_MMPBSA_sh_name

def mmpbsa_in()->list:
    mmpbsa_in_data = [
f"""
sample input file for running alanine scanning
 &general
   startframe=1, endframe=2000, interval=5, 
   verbose=1
/
&gb
  igb=66, saltcon=0.145
/
&pb
    istrng=0.145,
/
&alanine_scanning
/
"""] 
    return mmpbsa_in_data

def mmbpsa_in_gen(input_dict ) : 
    if input_dict["MMPBSA.IN_PATH"] == []:
        mmpbsa_in_list = mmpbsa_in()
        with open(f"mmpbsa.in", "w+") as mmpbsa:
            for line in mmpbsa_in_list :
                mmpbsa.write(f"{line}")
            mmpbsa.close()
            

def all_process_sh_gen(pdbfh_base_name, amber_source : str = "amber22/amber.sh"):
    all_process_sh = [f'''#!/bin/bash
#SBATCH --job-name=MD_pipeline
#SBATCH --partition=gpu
#SBATCH --exclusive
#SBATCH --output=MD_pipeline.out
#SBATCH --error=MD_pipeline.error
#SBATCH --time=72:00:00

echo "Loading modules..."
module load amber  
source {amber_source}

tleap -s -f tleap_mut.in > tut_leap.out

$AMBERHOME/bin/pmemd.cuda -O -i min.in -o min.out -p {pdbfh_base_name}_solvated.prmtop -c {pdbfh_base_name}_solvated.inpcrd \
-r min.rst -ref {pdbfh_base_name}_solvated.inpcrd
echo "min done"
$AMBERHOME/bin/pmemd.cuda -O -i heat.in -o heat.out -p {pdbfh_base_name}_solvated.prmtop -c min.rst \
-r heat.rst -x heat.mdcrd -ref min.rst

echo "heat done"
$AMBERHOME/bin/pmemd.cuda -O -i density.in -o density.out -p {pdbfh_base_name}_solvated.prmtop -c heat.rst \
-r density.rst -x density.mdcrd -ref heat.rst

echo "density done"
$AMBERHOME/bin/pmemd.cuda -O -i equil.in -o equil.out -p {pdbfh_base_name}_solvated.prmtop -c density.rst \
-r equil.rst -x equil.mdcrd

mpirun -v -np 4 $AMBERHOME/bin/pmemd.cuda.MPI -O -i prod.in -o prod1.out \
-p {pdbfh_base_name}_solvated.prmtop -c equil.rst -r prod1.rst -x prod1.mdcrd

echo "prod 1 done ..."
mpirun -v -np 4 $AMBERHOME/bin/pmemd.cuda.MPI -O -i prod.in -o prod2.out \
-p {pdbfh_base_name}_solvated.prmtop -c prod1.rst -r prod2.rst -x prod2.mdcrd

echo "prod 2 done..."
mpirun -v -np 4 $AMBERHOME/bin/pmemd.cuda.MPI -O -i prod.in -o prod3.out \
-p {pdbfh_base_name}_solvated.prmtop -c prod2.rst -r prod3.rst -x prod3.mdcrd

echo "prod 3 done..."
mpirun -v -np 4 $AMBERHOME/bin/pmemd.cuda.MPI -O -i prod.in -o prod4.out \
-p {pdbfh_base_name}_solvated.prmtop -c prod3.rst -r prod4.rst -x prod4.mdcrd''']
    
    
    
    with open("all_process.sh", "w+") as all_process_sh_file : 
        for line in all_process_sh : 
            all_process_sh_file.write(f"{line}\n")
        all_process_sh_file.close()
        all_process_sh_name = "all_process.sh"
        

def MPBSA_sh_helper( pdbfh_base_name, amber_source : str = "amber22/amber.sh") :
    '''
    pdbfh_base_name     : base name of the pdb file
    file_handle_mut_all : base bame of the mutated file 

    returns             : .sh file 
    '''
    
    mut_bash_sh = [f"#!/bin/bash",
        f"#SBATCH --job-name=run_66_mut",
        f"#SBATCH --partition=cpu",
        f"#SBATCH --output=run_mmpbsa_66.out",
        f"#SBATCH --error=run_mmpbsa_66.error",
        f"#SBATCH --time=48:00:00",
        f'''echo "Loading modules..."'''  ,  
        f"module load amber " ,
        f"source {amber_source}",
        f"""$AMBERHOME/bin/MMPBSA.py -O -i \
mmpbsa.in -o \
FINAL_RESULTS_MMPBSA_{pdbfh_base_name}.dat\
 -sp {pdbfh_base_name}_solvated.prmtop\
 -cp {pdbfh_base_name}.prmtop\
 -rp {pdbfh_base_name}_recpt.prmtop\
 -lp {pdbfh_base_name}_ligand.prmtop\
 -y *.mdcrd"""
        ]
    with open("run_MMPBSA.sh", "w+") as mut_bash_sh_file : 
            for line in mut_bash_sh : 
                mut_bash_sh_file.write(f"{line}\n")
            mut_bash_sh_file.close()
            mut_bash_sh_file_name = "run_MMPBSA.sh"

