#!/usr/bin/env python

from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import Lipinski
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Crippen

# Read in the reformatted sdf file directly using RDKit.
suppl = Chem.SDMolSupplier("./Sample_Out/exp_sdf_reformatted.sdf")

# Convert each object in the supplier to mol objects. tqdm is a module that
# shows a progress bar when wrapped around an iterable object.
print("Converting sdf compounds to mol objects.")
mols = [x for x in tqdm(suppl)]

# Sets up a list to hold the index of compounds that violate two of Lipinski's
# rules of five. The counter is to track that index number.
lipinski_violators = []
counter = 0

print("Scanning molecules for Lipinski violations.")
for mol in tqdm(mols):
    # Assume no violations
    dono_viol = False
    acceptor_viol = False
    mw_viol = False
    logp_viol = False

    # Use RDKit functions to get hdonors, acceptors, molecular weight and
    # logP.
    hdonors = Lipinski.NHOHCount(mol)
    hacceptors = Lipinski.NOCount(mol)
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    logp = Crippen.MolLogP(mol)

    # Make the checks if the current mol actually violates a role.
    if hdonors > 5:
        dono_viol = True
    if hacceptors > 10:
        acceptor_viol = True
    if mw > 500:
        mw_viol = True
    if logp > 5:
        logp_viol = True

    # Check if the violation sum is greater than one and assign the molecule
    # as a violator saving the index to a list.
    if sum([dono_viol, acceptor_viol, mw_viol, logp_viol]) > 1:
        lipinski_violators.append(counter)

    counter += 1

# Print the number of lipinski violators.
print(f"{len(lipinski_violators)} compounds violate more than one of", 
        "Lipinski's rules.")

# Write out these indices for later processing and removal from the sdf
# This script takes a long time to run so it's best to split up the script.
if len(lipinski_violators) > 0:
    print("Saving the bad molecule locations.")
    with open("./Sample_Out/strict_bad_mol_indices.txt", 'w') as f:
        for idx in lipinski_violators:
            f.write(f"{idx}\n")
