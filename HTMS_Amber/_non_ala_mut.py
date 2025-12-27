import os
import re

import _ala_mut
import _defaults
import _utils
from modeller import *
from modeller.automodel import autosched
from modeller.optimizers import ConjugateGradients, MolecularDynamics

""" The bullk of the code is taken from the Modeller example script
mutate_model.py, which is part of the Modeller distribution."""

def optimize(atmsel, sched):
    #conjugate gradient
    for step in sched:
        step.optimize(atmsel, max_iterations=200, min_atom_shift=0.001)
    #md
    refine(atmsel)
    cg = ConjugateGradients()
    cg.optimize(atmsel, max_iterations=200, min_atom_shift=0.001)
#molecular dynamics
def refine(atmsel):
    # at T=1000, max_atom_shift for 4fs is cca 0.15 A.
    md = MolecularDynamics(cap_atom_shift=0.39, md_time_step=4.0,
                           md_return='FINAL')
    init_vel = True
    for (its, equil, temps) in ((200, 20, (150.0, 250.0, 400.0, 700.0, 1000.0)),
                                (200, 600,
                                 (1000.0, 800.0, 600.0, 500.0, 400.0, 300.0))):
        for temp in temps:
            md.optimize(atmsel, init_velocities=init_vel, temperature=temp,
                         max_iterations=its, equilibrate=equil)
            init_vel = False

#use homologs and dihedral library for dihedral angle restraints
def make_restraints(mdl1, aln):
   rsr = mdl1.restraints
   rsr.clear()
   s = Selection(mdl1)
   for typ in ('stereo', 'phi-psi_binormal'):
       rsr.make(s, restraint_type=typ, aln=aln, spline_on_site=True)
   for typ in ('omega', 'chi1', 'chi2', 'chi3', 'chi4'):
       rsr.make(s, restraint_type=typ+'_dihedral', spline_range=4.0,
                spline_dx=0.3, spline_min_points = 5, aln=aln,
                spline_on_site=True)
       
#first argument
def _extract_mut_info(s):
    pattern = r'([A-Za-z])(\d+)([A-Za-z])'
    match = re.search(pattern, s)
    if match:
        respos = int(match.group(2))
        restyp_single_letter = str(match.group(3))
        restyp = _utils.amino_acids[restyp_single_letter]
    else:
        raise ValueError(f"Could not extract residue number from string: {s}")
    return respos, restyp

class mutate_model:
    def __init__(self, modelname, mutation , chain):
        self.modelname = modelname # #input file name, base pdb
        self.mutation = mutation
        self.chain = chain
        self.respos, self.restyp=_extract_mut_info(self.mutation)
        
        
    def _do_mut(self):
        self.s = Selection(self.mdl1.chains[self.chain].residues[self.respos])
        self.s.mutate(residue_type=self.restyp)
        
        
    def mutate_pdb(self):
       
        # Set a different value for rand_seed to get a different final model
        self.env = Environ(rand_seed=-49837)
        self.env.io.hetatm = True
        #soft sphere potential
        self.env.edat.dynamic_sphere=False
        #lennard-jones potential (more accurate)
        self.env.edat.dynamic_lennard=True
        self.env.edat.contact_shell = 4.0
        self.env.edat.update_dynamic = 0.39
        # Read customized topology file with phosphoserines (or standard one)
        self.env.libs.topology.read(file='$(LIB)/top_heav.lib')
        # Read customized CHARMM parameter library with phosphoserines (or standard one)
        self.env.libs.parameters.read(file='$(LIB)/par.lib')
        # Read the original PDB file and copy its sequence to the alignment array:
        self.mdl1 = Model(self.env, file=self.modelname)
        self.ali = Alignment(self.env)
        self.ali.append_model(self.mdl1,
                              atom_files=self.modelname, 
                              align_codes=self.modelname)

        self._do_mut()
        
        self.ali.append_model(self.mdl1, align_codes=self.modelname)
        self.mdl1.clear_topology()
        self.mdl1.generate_topology(self.ali[-1])
        self.mdl1.transfer_xyz(self.ali)
        self.mdl1.build(initialize_xyz=False, build_method='INTERNAL_COORDINATES')
        self.mdl2 = Model(self.env, file=self.modelname)
        self.mdl1.res_num_from(self.mdl2,self.ali)
        self.mdl1.write(file=self.modelname+self.mutation+'.tmp')
        self.mdl1.read(file=self.modelname+self.mutation+'.tmp')
        make_restraints(self.mdl1, self.ali)
        self.mdl1.env.edat.nonbonded_sel_atoms=1
        self.sched = autosched.loop.make_for_model(self.mdl1)
        
        
        self.s = Selection(self.mdl1.chains[self.chain].residues[self.respos])
        self.mdl1.restraints.unpick_all()
        self.mdl1.restraints.pick(self.s)
        self.s.energy()
        self.s.randomize_xyz(deviation=4.0)
        self.mdl1.env.edat.nonbonded_sel_atoms=2
        optimize(self.s, self.sched)
        self.mdl1.env.edat.nonbonded_sel_atoms=1
        optimize(self.s, self.sched)
        self.s.energy()
        self.mdl1.write(file=self.modelname+self.mutation+'.pdb')
        os.remove(self.modelname+self.mutation+'.tmp')
        
class multi_mutate_model(mutate_model):
    def __init__(self, modelname, mutation,  chain):
        self.modelname = modelname #input file name, base pdb
        self.chain = chain
        self.mutation = mutation 
        self.new_tuple_array_inputs = []
        for muts in self.mutation.split("_"):

            self.respos , self.restyp = _extract_mut_info(muts)
            self.new_tuple_array_inputs.append((self.respos, self.restyp))
        
    def _do_mut(self):
        for respos, restyp in self.new_tuple_array_inputs:
            s = Selection(self.mdl1.chains[self.chain].residues[respos])
            s.mutate(residue_type=restyp)
        
      
def general_mutate(modelname : str,
                   mutation : str,
                 chain = "A" ):
    
    if len(modelname.split("_")) >1  : 
        mutator= multi_mutate_model(modelname,mutation, chain)
        mutator.mutate_pdb()
    elif _extract_mut_info(modelname) is not None:
        mutator = mutate_model(modelname,mutation, chain)
        mutator.mutate_pdb()
    else :
        raise ValueError(f"Invalid modelname format: {modelname}. Expected format is '<resnum><resname>' or '<resnum><resname>_<resnum><resname>...'")
    

def split(pdbfh, pdbfh_base_name) :
    with open(pdbfh, "r") as f :
        pdb_data = f.readlines()
    #splits
    struct_pdb_data = _ala_mut.pdb_split(pdb_data, 0 )
    file_handle_structure = pdbfh_base_name + "_recpt.pdb"
    with open(file_handle_structure, "w+") as pdb_file : 
        for line in struct_pdb_data : 
            pdb_file.write(f"{line}")
        pdb_file.close()
    ligand_pdb = _ala_mut.pdb_split(pdb_data, 1 )
    file_handle_ligand = pdbfh_base_name + "_ligand.pdb"
    with open(file_handle_ligand, "w+") as pdb_file : 
        for line in ligand_pdb : 
            pdb_file.write(f"{line}")
        pdb_file.close()

def tleap_in_gen( pdbfh_base_name ): 
   
	tleap_mut_in = tleap_gen(pdbfh_base_name )
	with open("tleap_mut.in", "w+") as tleap : 
		for line in tleap_mut_in : 
			tleap.write(f"{line}\n")
		tleap.close()
	os.system(f"dos2unix tleap_mut.in") #not sure if needed. 
	tleap_file_name ="tleap_mut.in"
    
	return  tleap_file_name

def tleap_gen(pdbfh_base_name ) -> list:
    '''
    pdbfh_base_name     : base name of the pdb file
    returns             : tleap file as a list 
 
    
    '''

    tleap_in = [f"source leaprc.protein.ff19SB",
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
        f"quit"]
    return tleap_in


def non_ala_file_pop(pdbfh, amber_source ):
    pdbfh_base_name = os.path.basename(pdbfh).split(".")[0]
    split(pdbfh, pdbfh_base_name)
    tleap_file_name = tleap_in_gen(pdbfh_base_name)
    _defaults.all_process_sh_gen(pdbfh_base_name,amber_source = amber_source)
    _defaults.MPBSA_sh_helper(pdbfh_base_name, amber_source= amber_source)