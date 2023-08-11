
import pandas as pd
import os
import pybktree
import Levenshtein


######## Matching functions ########

# Exact matching
def exact_match(organism, substrate_df, termini_path, output_directory):
    termini = pd.read_csv(termini_path)
    substrate_subsites_filtered = substrate_df[substrate_df['organism'] == organism].reset_index(drop=True)
    for terminus in termini['Non-Tryptic Termini']:
        if len(terminus) != 8:
            raise ValueError(f"Invalid length of terminus: {terminus}. Each terminus should be of length 8.")
        
    matches = substrate_subsites_filtered[substrate_subsites_filtered['Cleavage Site'].isin(termini['Non-Tryptic Termini'])]
    matches = matches.drop(['Site_P4', 'Site_P3', 'Site_P2', 'Site_P1', 'Site_P1prime', 'Site_P2prime','Site_P3prime', 'Site_P4prime'], axis=1)
    
    temp_termini = termini.set_index('Non-Tryptic Termini') # set 'Non-Tryptic Termini' as index

    if len(matches) > 0:
        matches['Non-Tryptic Termini'] = matches['Cleavage Site']
        matches['Protein'] = matches['Non-Tryptic Termini'].map(temp_termini['Protein'])
        matches['Protein ID'] = matches['Non-Tryptic Termini'].map(temp_termini['Protein ID'])
        matches = matches.rename(columns={'cleavage_type': 'Cleavage Type', 'organism': 'Organism', 'Substrate_name': 'MEROPS Substrate'})
        matches = matches[['Non-Tryptic Termini', 'Protein', 'Protein ID', 'Cleavage Site', 'MEROPS Substrate', 'Protease', 'Cleavage Type',  'Organism']]
        # matches['Protein'] = matches['Cleavage Site'].map(termini.set_index('Non-Tryptic Termini')['Protein'])
        # matches['Protein ID'] = matches['Cleavage Site'].map(termini.set_index('Non-Tryptic Termini')['Protein ID'])
    else:
        matches = pd.DataFrame(columns=['Non-Tryptic Termini', 'Protein', 'Protein ID', 'Cleavage Site', 'MEROPS Substrate', 'Protease', 'Cleavage Type',  'Organism'])
    
    # print(matches)
    matches.to_csv(f'{output_directory}/exact_matches.csv', index=False)
    return matches


# Levenshtein distance fuzzy-matching
def fuzzy_match(organism, substrate_df, termini_path, output_directory):
    termini = pd.read_csv(termini_path)
    substrate_subsites_filtered = substrate_df[substrate_df['organism'] == organism].reset_index(drop=True)
    input_termini = termini['Non-Tryptic Termini']
    substrates = substrate_subsites_filtered['Cleavage Site']
    termini = termini.reset_index(drop=True)
    matches = []

    for terminus in input_termini:
        if len(terminus) != 8:
            raise ValueError(f"This file contains termini with invalid lengths. Each terminus should be of length 8 to match to cleavage sites found in the MEROPS database.")
    # Define Levenshtein distance metric
    def metric(a, b): # a and b are pairs of index, string -> string is a non-tryptic terminus/cleavage site
        # Check exact match for the 4th and 5th amino acids
        if a[1][3:5] != b[1][3:5]: # if the 4th and 5th amino acids aren't the same, return a distance of inf (bad match)
            return float('inf')

        # Calculate Levenshtein distance for the surrounding parts (excluding 4th and 5th amino acids )
        surrounding_a = a[1][:3] + a[1][5:]
        surrounding_b = b[1][:3] + b[1][5:]
        return Levenshtein.distance(surrounding_a, surrounding_b) # return Levenshtein distance for surrounding amino acids

    tree = pybktree.BKTree(metric, list(enumerate(substrates)))

    for i, terminus in enumerate(input_termini):
        for _, (j, cleave_site) in tree.find((None, terminus), 2):  # Find matches -> set Levenshtein distance to 2
            matches.append((terminus, termini.loc[i,'Protein'], termini.loc[i,'Protein ID'], cleave_site, substrate_subsites_filtered.loc[j, 'Substrate_name'], substrate_subsites_filtered.loc[j,'Protease'], substrate_subsites_filtered.loc[j,'cleavage_type'],substrate_subsites_filtered.loc[j, 'organism']))

            
    matches_df = pd.DataFrame(matches, columns=['Non-Tryptic Termini', 'Protein', 'Protein ID', 'Cleavage Site', 'MEROPS Substrate', 'Protease', 'Cleavage Type',  'Organism'])
    matches_df = matches_df.drop_duplicates()
    matches_df.to_csv(f'{output_directory}/fuzzy_matches.csv', index=False)
    return matches_df

# organisms, substrate_subsites = substrate_processing('./substrate.csv')
# termini = pd.read_csv('./test2_output/non_tryp_termini.csv')
# fuzzy_match("Escherichia coli", substrate_subsites, termini, '.')
