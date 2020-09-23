#!/usr/bin/bash

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Read in processed csv here, from schrodinger_blind_screen_result_analysis.py 
# Data here from two different crystal structures.
hosu = pd.read_csv('new_chemotypes_RJ77_SP_results_processed.csv')
et = pd.read_csv('new_chemotypes_6ET4_nocofactors_SP_results_processed.csv')

# Rename columns to be congruent with one another.
hosu.columns = ['Compound', 'Docking Score']
et.columns = ['Compound', 'Docking Score']

# Tag compounds that have been suggested vs. those tested for later plotting.
# first is docking in crystal structure RJ77
hosu_category = []
for comp in hosu['Compound'].tolist():
    if 'AML' in comp:
        hosu_category.append('Newly Proposed')
    else:
        hosu_category.append('Already Tested')

# Repeat the same only for 6ET4
et_category = []
for comp in et['Compound'].tolist():
    if 'AML' in comp:
        et_category.append('Newly Proposed')
    else:
        et_category.append('Already Tested')

# Add the compound categories to the pandas dataframe
hosu['Compound Category'] = hosu_category
et['Compound Category'] = et_category

# Subset both dataframes showing only newly proposed for analysis
new_comps = et[et['Compound Category'] == 'Newly Proposed']

new_comps_hos = hosu[hosu['Compound Category'] == 'Newly Proposed']

# Print Results
print(new_comps)
print(new_comps_hos)

# See where lead target ends up
print(hosu[hosu['Compound'] == 'HOSU3'])
print(et[et['Compound'] == 'HOSU3'])

print(hosu.head(60))

# Plot results using R language ggplot style.
plt.style.use("ggplot")

# Set up seaborn swarmplot as an axis object, add labels, and save.
# This will only save structures docked in 6ET4.
ax = sns.swarmplot(y=et["Docking Score"], hue=et["Compound Category"],
        x=[""]*len(et),
        palette={"Already Tested": "steelblue",
            "Newly Proposed": "orange"})
ax.set_ylabel("Docking Score (kcal/mol)")
ax.set_title("Ranking of Proposed DHODH Inhibitors vs Tested Compounds" + 
        "\n" + "in 6ET4 Crystal")
plt.savefig("new_chemotypes_6et4_docked_scores.png", dpi=300, 
        bbox_inches="tight")
