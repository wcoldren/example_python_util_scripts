from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem

def make_smiles_dict(filename):
    """Takes in a csv file with two columns. The first column
    would be a compound identifier name while the second is a smiles
    string for the molecule

    This function returns a dictionary where the keys are compound
    names and values are the smiles stings."""

    smiles_raw = open(filename).read().splitlines()
    smiles_dict = {}

    for line in smiles_raw:
        smiles_dict[line.split(',')[0]] = line.split(',')[1]
    
    return smiles_dict

def smiles_dict_to_mol_list(smiles_dict):
    """smiles dict is a dictionary object containing molecule names
    as keys and smiles strings as values.

    The return value is a list of RDKit mol objects.
    """
    
    smiles_as_mol = []
    for mol_name, smiles_string in smiles_dict.items():
        try:
            mol = Chem.MolFromSmiles(smiles_string)
            mol.SetProp("_Name", mol_name)
            smiles_as_mol.append(mol)
        except:
            print("Error processing:", mol_name)

    return smiles_as_mol

def mol_to_sdf(mol_list: list, output: str):
    """Takes as input a list of mol objects
    as well as an output file name. The mol objects will be
    hydrogenated, embedded, and saved to an sdf format file."""

    hydrogenated = []
    for m in mol_list:
        m_h = Chem.AddHs(m)
        AllChem.EmbedMolecule(m_h)
        hydrogenated.append(m_h)
    # Create an SDF writer object and write all mol objects to an sdf file.
    w = Chem.SDWriter(output)
    for m in hydrogenated: w.write(m)

def fingerprint_mols(mols: list):
    """Accepts a list of mol objects and outputs a list of their RDKit
    fingerprints."""
    fps = [Chem.RDKFingerprint(x) for x in mols]

    return fps

if __name__ == '__main__':
    main()
