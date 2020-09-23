from docking_prep_tools import *
import matplotlib.pyplot as plt

# The three following lines utilize functions defined in the 
# docking_prep_tools.py script. Change variables as needed.
all_smiles = make_smiles_dict('sample_smiles.smi')

mols = smiles_dict_to_mol_list(all_smiles)

all_fps = fingerprint_mols(mols)

# Once fingerprinted the tanimoto comparison can be calculated and printed
# as follows, this is just comparing two molecules.
print(DataStructs.FingerprintSimilarity(all_fps[-1], all_fps[-2], 
    metric=DataStructs.TanimotoSimilarity))

# This loops through and compares all molecules each other.
# Adjust as needed, a slicker way to do this would be with itertools.
similarities = []
for i in range(len(all_fps)):
    for j in range(len(all_fps)):
        if j > i:
            similarities.append(DataStructs.FingerprintSimilarity(all_fps[i],
                all_fps[j]))

# Create a quick histogram plot to view results.
plt.style.use("ggplot")
plt.hist(similarities)
plt.xlim(0, 1)
plt.title("Tanimoto Similarity Indexes for 100 Sample Compounds")
plt.xlabel("Tanimoto Coefficient")
plt.ylabel("Number of Molecule Comparisons")

# Save figure as a png
plt.savefig("sample_tanimoto.png", dpi=300, bbox_inches="tight")
