from pubchempy import *
from rdkit import Chem
from rdkit.Chem import AllChem
import json, subprocess
import autode as ade
import pandas as pd

def get_LA_smiles(LA_file="LAs/raw_names.txt"):

    raw_names = []
    with open(LA_file, "r") as f:
        for line in f:
            raw_names.append(line.replace('\n',''))

    SMILES = []
    raw_ID_dict = {}
    LA_dict = {}
    for i, LA in enumerate(raw_names):
        c = get_compounds(LA, 'name')
        try:
            smiles = c[0].isomeric_smiles
            SMILES.append(smiles)
            raw_ID_dict[str(i)] = LA
            LA_dict[LA] = smiles
        except:
            print(f"error {LA}")

    with open("LAs/smiles.txt", "w") as f:
        for s in SMILES:
            f.write(s+"\n")

    LA_id_dump = open("LAs/LA_dict.json", "w")
    raw_id_dump = open("LAs/id.json", "w")
    json.dump(LA_dict, LA_id_dump)
    json.dump(raw_ID_dict, raw_id_dump)
    
    return raw_ID_dict


def make_complex(epoxide="C1C(O1)C2=CC=CC=C2", metal_symbols=['Al', 'Ag', 'B', 'Cu', 'Ga', 'Li', 'In'], oxygen_index=2):

    epoxide_mol = Chem.MolFromSmiles(epoxide)
    SMILES = []
    try:
        with open("LAs/smiles.txt", "r") as f:
            for line in f:
                line = line.replace('\n','')
                SMILES.append(line)
                ID_dict = json.load(open("LAs/id.json", "r"))
    except:
        SMILES, ID_dict = get_LA_smiles()

    complexes_ID_dict = {}
    list_of_complexes = []
    for i, s in enumerate(SMILES):
        LAMol = Chem.MolFromSmiles(s) # lewis acid rdkit mol object
        cplex = Chem.CombineMols(epoxide_mol, LAMol) # combine styrene oxide and lewis acid
        for atom in cplex.GetAtoms(): # iterate over atoms
            if atom.GetSymbol() in metal_symbols: # if atom is a metal
                metal_index = atom.GetIdx() # set metal index

        edcomplex = Chem.EditableMol(cplex) # editable molecule object
        edcomplex.AddBond(oxygen_index, metal_index, order=Chem.rdchem.BondType.SINGLE) # add bond
        cplex = edcomplex.GetMol() # editable molecule back to rdkit mol object

        smiles = Chem.MolToSmiles(cplex) # convert complex mol object back to smiles
        list_of_complexes.append(smiles)

        complexes_ID_dict[ID_dict[str(i)]] = smiles
        
    textfile = open("LAs/complexes.smi", "w")
    for element in list_of_complexes:
        textfile.write(element + "\n")
    textfile.close()

    id_dump = open("LAs/complex_id.json", "w")
    json.dump(complexes_ID_dict, id_dump)

def get_rmsd(ref, geom):

    cmd = ('obrms', f'{ref}', f'{geom}')
    rmsd_file = open("rmsd", "w")
    subprocess.call(cmd, stdout=rmsd_file)
    rmsd_file.close()
    with open("rmsd", "r") as f:
        lines = f.readlines()
        for line in lines:
            rmsd = line.split()[-1]

    return float(rmsd)
        
    
def benchmark_DFT(ID_dict={"name": None, "SMILES": None}, dft_dict=None, methods=None, solvent=None, df=None, n_cores=8):

    orca = ade.methods.ORCA() # set method
    solvent_block=f'%cpcm\nsmd true\nSMDsolvent "{solvent}"\nend' # solvent block
    energies = []
    ID = ID_dict["name"]
    SMILES = ID_dict["SMILES"]
    for method_name, keywords in dft_dict.items(): # iterate over defined methods
        if method_name in methods: # if this method is in the list of methods you want to calculate
            if "(" in method_name:
                method_name = method_name.replace("(", "")
                method_name = method_name.replace(")", "")
            if SMILES is not None:
                molecule = ade.Molecule(name=f'{ID}_{method_name}', smiles=SMILES)
            else:
                molecule = ade.Molecule(f"xtb_opts/{ID}_xtb_opt.xyz") # create autoDE molecule
            print(f"running {method_name}")
            optimization = ade.Calculation(name=f'{ID}_{method_name}',
                                            molecule=molecule,
                                            method=orca,
                                            keywords=keywords,
                                            n_cores=n_cores,
                                            other_input_block=solvent_block)
            optimization.generate_input()
            optimization.output.filename = f"{ID}_{method_name}.out"
            optimization.execute_calculation()
            optimization.clean_up()
            optimization._add_to_comp_methods()
            energy = optimization.get_energy()
            energies.append(energy.to("kcal mol-1"))
            print(f"done {method_name}")
            
    rmsd_dict = {}

    ref = f"{ID}_{methods[-1]}_orca.xyz"
    rmsd_dict[ref] = "ref"
    for method in methods[:-1]:
        geom = f"{ID}_{method}_orca.xyz"
        rmsd = get_rmsd(ref, geom)
        rmsd_dict[geom] = rmsd
    rmsd_lst = list(rmsd_dict.values())
    df[f"rmsd {ID}"] = rmsd_lst[::-1]
    energies = [i - energies[-1] for i in energies] 
    energies[-1] = "ref"
    df[f"energies {ID}"] = energies
    return df
    
        
        
