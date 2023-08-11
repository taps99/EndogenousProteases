import pandas as pd
import os
import numpy as np
import time
import re
import requests
from weblogo import *


peptide_df_all = pd.read_csv("./peptide_output_with_sequence.csv")


# subset of peptide_df that only contains peptides that appear in at least 3 of the 4 wild-type samples, but don't appear in mutant samples
subset_df_wt = peptide_df_all[(peptide_df_all.iloc[:, 7:11].sum(axis=1)  == 0) & (peptide_df_all.iloc[:, 11:].sum(axis=1) >= 3)]


# Now need to acquire protein sequence from uniprot using accession ID for each of these peptides
# get protein sequence based on UniProt accession ID
def retrieve_protein_sequence(id):
    # URL for FASTA file for a specific accession ID
    url = f"https://www.uniprot.org/uniprot/{id}.fasta" 
    response = requests.get(url)
    response.raise_for_status()
    sequence_lines = response.text.split('\n')[1:]  # Skip the header line starting with ">"
    sequence = ''.join(sequence_lines)
    # print(sequence)
    time.sleep(3)
    return sequence

# store sequence as new column in df
# subset_df_wt['Protein Sequence'] = subset_df_wt['Protein ID'].apply(retrieve_protein_sequence)
subset_df_wt.to_csv('./WT_only_subset_3.csv', index=False)

# subset_df_wt = pd.read_csv('./WT_only_subset_3.csv')
# Now for each non-tryptic cleavage, get the previous 5 amino acids before the cleavage, and the 5 amino acids following the cleavage
# maybe also make a version that shows where the cleavage actually occurs (with a '/' or '.') and put that in a column
# print(subset_df_wt_seq['Protein Sequence'].iterrows)


def peptide_seq_match(df):
    # match peptide with protein sequence using re 
    peptide = df['Peptide']
    prev_AA = df['Prev AA']
    tryp_state = df['Tryptic State']
    sequence = str(df['Protein Sequence'])
    matches = re.finditer(peptide, sequence)
    print(matches)

    n_terminal_list = []
    c_terminal_list = []

 
   
    for match in matches:
        # Gives start and end indices that you can use to get the amino acids before and after the match in a given sequence
        start_index = match.start()
        end_index = match.end()

        # Calculate indices for extracting the 5 amino acids before and after the match
        # 2nd argument acts as a boundary, so there's no values that go outside of the sequence length that can't be indexed (can't go below 0, can't go above length of sequence)
        before_start_index = max(start_index - 5, 0)
        after_end_index = min(end_index + 5, len(sequence))

        # Extract 5 amino acids before and after the match
        before_match = sequence[before_start_index:start_index]
        after_match = sequence[end_index:after_end_index]

        n_terminal = before_match + peptide[:5]
        c_terminal = peptide[-5:] + after_match
        
        if tryp_state == 'Non-Tryptic':
            n_terminal_list.append(n_terminal)
            c_terminal_list.append(c_terminal)
        else:
            if prev_AA == 'K' or prev_AA == 'R':
                c_terminal_list.append(c_terminal)
            elif peptide.endswith('K') or peptide.endswith('R'):
                n_terminal_list.append(n_terminal)
    # Add dashes to the end of sequences that are less than 10 aa long for N-terminal cleavages
    for i in range(len(n_terminal_list)):
        string = n_terminal_list[i]
        while len(string) < 8:
            string = "-" + string
            
        n_terminal_list[i] = string
    # Add dashes to beginning of sequences that are less than 10 aa long for C-terminal cleavages
    for i in range(len(c_terminal_list)):
        string = c_terminal_list[i]
        while len(string) < 8:
            string += "-"
        c_terminal_list[i] = string
            
    n_terminal_str = ''.join(n_terminal_list)
    c_terminal_str = ''.join(c_terminal_list)
    
    terminal_series = pd.Series({'N-terminal cleavage': n_terminal_str, 'C-terminal cleavage': c_terminal_str})
    terminal_series.to_csv('blahblah.csv', index=False)
    return terminal_series


subset_df_wt[['N-terminal cleavage', 'C-terminal cleavage']] = subset_df_wt.apply(peptide_seq_match, axis=1)

subset_df_wt.to_csv('./WT_only_subset_3.csv', index=False)

# Alignment?

# subset_df_wt['C-terminal cleavage'].to_csv('testo.csv', index=False, header=False)
subset_df_wt[subset_df_wt['C-terminal cleavage'] != '']['C-terminal cleavage'].to_csv('c_terminal.csv', index=False)
subset_df_wt[subset_df_wt['N-terminal cleavage'] != '']['N-terminal cleavage'].to_csv('n_terminal.csv', index=False)
c_terminal_df = subset_df_wt[subset_df_wt['C-terminal cleavage'] != '']['C-terminal cleavage']
n_terminal_df = subset_df_wt[subset_df_wt['N-terminal cleavage'] != '']['N-terminal cleavage']

non_tryp_seqs = pd.concat([c_terminal_df, n_terminal_df])
non_tryp_seqs = pd.DataFrame(non_tryp_seqs, columns=['Non-Tryptic Sites'])

# print((non_tryp_seqs))
# for index, row in subset_df_wt.iterrows():
#     if row['C-terminal cleavage'] in non_tryp_seqs['Non-Tryptic Sites']:
#         print(row['C-terminal cleavage'])
#         non_tryp_seqs['Protein ID'] = subset_df_wt['Protein ID']

# Put all non-tryptic C-termini and N-termini in one column
for index, row in non_tryp_seqs.iterrows():
    sites = row['Non-Tryptic Sites']
    if sites in subset_df_wt['C-terminal cleavage'].values:
        protein_id = subset_df_wt.loc[subset_df_wt['C-terminal cleavage'] == sites, 'Protein ID'].iloc[0]
        non_tryp_seqs.at[index, 'Protein ID'] = protein_id
    elif sites in subset_df_wt['N-terminal cleavage'].values:
        protein_id = subset_df_wt.loc[subset_df_wt['N-terminal cleavage'] == sites, 'Protein ID'].iloc[0]
        non_tryp_seqs.at[index, 'Protein ID'] = protein_id


non_tryp_seqs.to_csv('non_tryptic_sequences.csv', index=False)
# Sequence logo to visualize
# read_seq