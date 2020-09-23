#!/usr/bin/env python
import argparse
import pandas as pd
import numpy as np

def get_glide_docking_score(glide_csv_file, remove=True):
    '''This function will test a glide csv file and parse ligand and docking
    scores (this includes epik penalization modified from the basic Glide
    score. Duplicate ligands will be removed automatically unless flag turned
    off.'''
    
    # Check to make sure file passed is actually a schrodinger csv.
    if not glide_csv_file.endswith('.csv'):
        return glide_csv_file + 'is not a csv file!'
        
    my_data = pd.read_csv(glide_csv_file, low_memory=False)
    if np.isnan(my_data.loc[0, 'docking score']):
        name_series = my_data.loc[1:, 'Title']
        score_series = my_data.loc[1:, 'docking score']
        score_df = pd.concat([name_series, score_series], axis=1)
    else:
        name_series = my_data['Title']
        score_series = my_data['docking score']
        score_df = pd.concat([name_series, score_series], axis=1)
    
    # Check to see if there is an error in number of scores vs. number of
    # ligands.
    if len(score_df['Title'].tolist()) != len(score_df['docking score']):
        return glide_csv_file + 'has a mismatch in structures and scores!'
    
    # Remove duplicate ligands and reindex starting at 0
    if remove:
        score_df.drop_duplicates(subset='Title', inplace=True)
        score_df.reset_index(inplace=True)
    
    return score_df.loc[:, 'Title':]

def parse_arguments():
    parser = argparse.ArgumentParser(description=
            'Read in a blind screening Glide csv file for parsing.')
    parser.add_argument('--csv', type=str, nargs='+', help='Input Glide csv file.')

    args = parser.parse_args()
    return args


def main(input_files):
    
    for csv in input_files:
        file_base = csv.split('.')[0]
        result_df = (get_glide_docking_score(csv))
   
        result_df.to_csv('%s_processed.csv' % file_base, index=False)
   
if __name__ == '__main__':
    arguments = parse_arguments()
    main(arguments.csv)
