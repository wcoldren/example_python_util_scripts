#!/usr/bin/env python
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem.FilterCatalog import *
import os.path

if os.path.isfile("./Input/strict_bad_mol_indices.txt"):
    lipinski = open("./Input/strict_bad_mol_indices.txt").read().splitlines()
    # Overwrite the lipinski list with a deep copy of integer values for the
    # indices read in (they default as strings).
    lipinski[:] = [int(x) for x in lipinski]

    suppl = Chem.SDMolSupplier("./Input/exp_sdf_reformatted.sdf")
    print("converting chembridge SDF to mol objects.")
    mols = [x for x in tqdm(suppl)]

    print("Removing Lipinski violations...")
    lip_removed = []
    for i in range(len(mols)):
        if i not in lipinski:
            lip_removed.append(mols[i])

else:
    lip_removed = Chem.SDMolSupplier("./Input/exp_sdf_reformatted.sdf")
    print("There were no bad molecular indices, continuing to PAINS")

# Set up the RDKit PAINS filter
params = FilterCatalogParams()
params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
catalog = FilterCatalog(params)

print("Locating and removing PAINS")
# Set up a counter to track number of pains compounds and an empty list
# to store the good molecules.
pains = 0
good_molecules = []

# The best way to identify which molecules are pains is to use try and
# accept syntax. You try to get the PAINS first match saved as an entry, 
# but if the molecule is not pains entry becomes a nonetype object. If you
# run GetProp on this the program will exit on error. So if they molecule is
# good you hit the except clause and append the molecule to the good molecules
# list.
for i in tqdm(range(len(lip_removed))):
    entry = catalog.GetFirstMatch(lip_removed[i])
    try:
        entry.GetProp("Scope")
        pains += 1
    except:
        good_molecules.append(lip_removed[i])

# Write out the PAINS free and Lipinski filtered moleclues to an sdf file.
print("Writing output")
w = Chem.SDWriter("./Sample_Out/clean_exp_mols.sdf")
for m in good_molecules: w.write(m)

