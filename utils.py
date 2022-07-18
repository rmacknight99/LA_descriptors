from pubchempy import *
from rdkit import Chem
from rdkit.Chem import AllChem
import json

def get_LA_smiles(LA_file="LAs/raw_names.txt"):

    raw_names = []
    with open(LA_file, "r") as f:
        for line in f:
            raw_names.append(line.replace('\n',''))

    SMILES = []
    raw_ID_dict = {}
    for i, LA in enumerate(raw_names):
        c = get_compounds(LA, 'name')
        try:
            smiles = c[0].isomeric_smiles
            SMILES.append(smiles)
            raw_ID_dict[str(i)] = LA
        except:
            print(f"error {LA}")

    with open("LAs/smiles.txt", "w") as f:
        for s in SMILES:
            f.write(s+"\n")

    id_dump = open("LAs/id.json", "w")
    json.dump(raw_ID_dict, id_dump)
    
    return SMILES, raw_ID_dict

            
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

        complexes_ID_dict[smiles] = ID_dict[str(i)]

    textfile = open("LAs/complexes.smi", "w")
    for element in list_of_complexes:
        textfile.write(element + "\n")
    textfile.close()

    id_dump = open("LAs/complex_id.json", "w")
    json.dump(complexes_ID_dict, id_dump)

make_complex()
