
import pandas as pd
import os
import pybktree
import Levenshtein


    # peptide_df = pd.read_csv(df_path)
    # peptide_df = pd.read_csv('Good_Output/nontryp_pept_sequences.csv')
    # peptide_df = pd.read_csv('test2_output/non-tryp_pept_sequences.csv')

    # nontryp_termini = pd.concat([peptide_df['N-terminal'], peptide_df['C-terminal']])
    # nontryp_termini = pd.DataFrame(nontryp_termini, columns=['Non-Tryptic Termini'])
    # nontryp_termini_only = nontryp_termini.dropna()

    # for i, row in nontryp_termini.iterrows():
    #     sites = row['Non-Tryptic Termini']
    #     if sites in peptide_df['C-terminal'].values:
    #         nontryp_termini.at[i, 'Protein ID'] = peptide_df.loc[peptide_df['C-terminal'] == sites, 'Protein ID'].iloc[0]
    #         nontryp_termini.at[i, 'Protein'] = peptide_df.loc[peptide_df['C-terminal'] == sites, 'Protein'].iloc[0]
    #     elif sites in peptide_df['N-terminal'].values:
    #         nontryp_termini.at[i, 'Protein ID'] = peptide_df.loc[peptide_df['N-terminal'] == sites, 'Protein ID'].iloc[0]
    #         nontryp_termini.at[i, 'Protein'] = peptide_df.loc[peptide_df['N-terminal'] == sites, 'Protein'].iloc[0]
            

    # nontryp_termini = nontryp_termini.dropna()
    # nontryp_termini.to_csv('non_tryp_termini.csv', index=False)

def substrate_processing(path):

    ### Substrate table processing
    amino_acid_dict = {
        'Ala': 'A',
        'Arg': 'R',
        'Asn': 'N',
        'Asp': 'D',
        'Cys': 'C',
        'Gln': 'Q',
        'Glu': 'E',
        'Gly': 'G',
        'His': 'H',
        'Ile': 'I',
        'Leu': 'L',
        'Lys': 'K',
        'Met': 'M',
        'Phe': 'F',
        'Pro': 'P',
        'Ser': 'S',
        'Thr': 'T',
        'Trp': 'W',
        'Tyr': 'Y',
        'Val': 'V'
    }

    # current_dir = os.path.dirname(os.path.realpath(__file__))
    # csv_path = os.path.join(current_dir, 'substrate.csv')
    substrate_df = pd.read_csv(path)

    # Subset substrate_df to include relevant columns
    substrate_subsites = substrate_df.loc[:,['Substrate_name', 'Site_P4', 'Site_P3', 'Site_P2', 'Site_P1', 'Site_P1prime', 'Site_P2prime','Site_P3prime', 'Site_P4prime', 'organism', 'Protease', 'cleavage_type']]

    # Replace abbreviations with single letter code for amino acids
    substrate_subsites = substrate_subsites.replace(amino_acid_dict)

    # Join the amino acids in each subsite together to create an amino acid sequence representing the cleavage site of a substrate 
    substrate_subsites['Cleavage Site'] = substrate_subsites[['Site_P4', 'Site_P3', 'Site_P2', 'Site_P1', 'Site_P1prime', 'Site_P2prime','Site_P3prime', 'Site_P4prime']].apply(lambda row: ''.join(row.values.astype(str)), axis=1) # 95653 substrates total
    substrate_subsites = substrate_subsites.dropna(subset=['organism', 'Substrate_name', 'Protease'])
    # print(len(substrate_subsites)) # 87968 substrates after dropping NA's

    # Organism list
    organisms = sorted(substrate_subsites['organism'].str.lstrip().str.capitalize().unique().tolist()) # 1417 unique organisms -> 1397 unique organisms
    org = pd.Series(organisms)
    # org.to_csv('organisms.csv', index=False, header=False)
    return organisms, substrate_subsites

# Exact matching
def exact_match(organism, substrate_df, termini, output_directory):
    substrate_subsites_filtered = substrate_df[substrate_df['organism'] == organism].reset_index(drop=True)
    matches = substrate_subsites_filtered[substrate_subsites_filtered['Cleavage Site'].isin(termini['Non-Tryptic Termini'])]
    matches = matches.drop(['Site_P4', 'Site_P3', 'Site_P2', 'Site_P1', 'Site_P1prime', 'Site_P2prime','Site_P3prime', 'Site_P4prime'], axis=1)
    # print(matches)
    temp_termini = termini.set_index('Non-Tryptic Termini') # set 'Non-Tryptic Termini' as index
    # print(temp_termini)  
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

# organisms, substrate_subsites = substrate_processing('./substrate.csv')
# termini = pd.read_csv('./Good_Output/non-tryp_termini.csv')
# exact_match('Klebsiella pneumoniae', substrate_subsites,termini,'.')




# Levenshtein distance fuzzy-matching
def fuzzy_match(organism, substrate_df, termini, output_directory):
    substrate_subsites_filtered = substrate_df[substrate_df['organism'] == organism].reset_index(drop=True)
    input_termini = termini['Non-Tryptic Termini']
    substrates = substrate_subsites_filtered['Cleavage Site']
    termini = termini.reset_index(drop=True)
    matches = []

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
