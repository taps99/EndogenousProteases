import pandas as pd
import os
import numpy as np
import time
import re
import requests
from weblogo import *


peptide_df_all = pd.read_csv("./peptide_output.csv")


# subset of peptide_df that only contains peptides that appear in at least 3 of the 4 wild-type samples, but don't appear in mutant samples
subset_df_wt = peptide_df_all[(peptide_df_all.iloc[:, 7:11].sum(axis=1)  == 0) & (peptide_df_all.iloc[:, 11:].sum(axis=1) >= 2)]


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
# subset_df_wt.to_csv('./WT_only_subset_2.csv')

subset_df_wt = pd.read_csv('./WT_only_subset_2.csv')
# Now for each non-tryptic cleavage, get the previous 5 amino acids before the cleavage, and the 5 amino acids following the cleavage
# maybe also make a version that shows where the cleavage actually occurs (with a '/' or '.') and put that in a column
# print(subset_df_wt_seq['Protein Sequence'].iterrows)


def peptide_seq_match(df):
    # match peptide with protein sequence using re 
    peptide = df['Peptide']
    sequence = str(df['Protein Sequence'])
    matches = re.finditer(peptide, sequence)

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
        
        n_terminal_list.append(n_terminal)
        c_terminal_list.append(c_terminal)


    # Add dashes to the end of sequences that are less than 10 aa long for N-terminal cleavages
    for i in range(len(n_terminal_list)):
        string = n_terminal_list[i]
        while len(string) < 10:
            string = "-" + string
            
        n_terminal_list[i] = string
    # Add dashes to beginning of sequences that are less than 10 aa long for C-terminal cleavages
    for i in range(len(c_terminal_list)):
        string = c_terminal_list[i]
        while len(string) < 10:
            string += "-"
        c_terminal_list[i] = string
            
    n_terminal_str = ''.join(n_terminal_list)
    c_terminal_str = ''.join(c_terminal_list)

    return pd.Series({'N-terminal cleavage': n_terminal_str, 'C-terminal cleavage': c_terminal_str})


subset_df_wt[['N-terminal cleavage', 'C-terminal cleavage']] = subset_df_wt.apply(peptide_seq_match, axis=1)

subset_df_wt.to_csv('./WT_only_subset_2.csv', index=False)

# Alignment?

subset_df_wt['C-terminal cleavage'].to_csv('testo.csv', index=False, header=False)

# Sequence logo to visualize
# read_seq