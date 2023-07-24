import pandas as pd
import os
import requests
import concurrent.futures
import re
import numpy as np
import threading


def process_data(listbox_entries):
    unprocessed_dfs = []
    processed_dfs = []
    count = 0

    for entry in listbox_entries:
        df = pd.read_csv(entry, delimiter='\t', low_memory=False)
        unprocessed_dfs.append(df)

        condition1 = ~(df['Peptide'].str.endswith('K') | df['Peptide'].str.endswith('R') & df['Prev AA'].isin(['K', 'R']))
        condition2 = ~(df['Prev AA'].isin(['K', 'R']) & df['Next AA'].isin(['-']))
        condition3 = ~(df['Prev AA'].isin(['-']) & (df['Peptide'].str.endswith('K') | df['Peptide'].str.endswith('R')))

        merged_df = df[condition1 & condition2 & condition3]
        merged_df = merged_df.drop_duplicates().reset_index(drop=True)

        prev_aa_mask = merged_df['Prev AA'].str.contains('K|R')
        peptide_end_mask = merged_df['Peptide'].str.endswith('K') | merged_df['Peptide'].str.endswith('R')

        merged_df.loc[prev_aa_mask | peptide_end_mask, 'Tryptic State'] = 'Semi-Tryptic'
        merged_df.loc[~(prev_aa_mask | peptide_end_mask), 'Tryptic State'] = 'Non-Tryptic'

        merged_df["Peptide:Protein"] = merged_df[['Peptide', 'Protein Description']].agg(':'.join, axis=1)
        merged_df_subset = merged_df[['Peptide:Protein', 'Tryptic State', 'Peptide', 'Prev AA', 'Next AA', 'Protein Description', 'Protein ID']]
        
        processed_dfs.append(merged_df_subset) # List of dataframes containing non-tryptic peptides, where a dataframe represents one sample
    return processed_dfs, unprocessed_dfs


def output_files(processed_filepaths, processed_dfs, unprocessed_dfs, output_directory):

    processed_dict = dict(zip(processed_filepaths, processed_dfs))
    unprocessed_dict = dict(zip(processed_filepaths, unprocessed_dfs))


    # This is to get the file that contains ALL peptides from the input files (includes tryptic and non-tryptic)
    combined_df = pd.concat(unprocessed_dfs, ignore_index=True)
    combined_df["Peptide:Protein"] = combined_df[['Peptide', 'Protein Description']].agg(':'.join, axis=1)
    combined_df.to_csv(os.path.join(output_directory, 'all_peptides.csv'), index=False)


    # This is to get the unique proteins from the entire dataset (includes tryptic and non-tryptic)
    combined_df_subset = combined_df[['Peptide:Protein', 'Peptide', 'Prev AA', 'Next AA', 'Protein Description', 'Protein ID']].drop_duplicates(subset='Protein ID')
    
    for i, df in enumerate(unprocessed_dict.values()):
        sample_name = list(unprocessed_dict.keys())[i]
        # if len(sample_name) > 2:
        for protein in combined_df_subset:
            combined_df_subset[sample_name] = combined_df_subset['Protein ID'].isin(df['Protein ID']).astype(int)

    combined_df_subset.to_csv(os.path.join(output_directory, 'all_unique_proteins.csv'), index=False)

    proteases = combined_df_subset[combined_df_subset['Protein Description'].str.contains(r'(?:protease|peptidase)')]
    proteases.to_csv(os.path.join(output_directory, 'all_proteases.csv'), index=False)
           

    # This part involves the actual output of non-tryptic peptides + relevant information
    unique_peptide_dfs=[]
    for df in processed_dfs: # GETS UNIQUE PEPTIDES FROM EACH SAMPLE
        unique_peptides = df.drop_duplicates(subset=['Peptide:Protein'])
        unique_peptide_dfs.append(unique_peptides) 

    combined_peptide_df = pd.concat(unique_peptide_dfs) # CONCAT THE NON-TRYPTIC PEPTIDES FROM EACH SAMPLE INTO ONE LARGE DF
    # THEN ONLY KEEP UNIQUE ONES, WHICH MEANS THIS ONE CONTAINS ALL UNIQUE PEPTIDES FROM ALL SAMPLES
    # print(len(combined_peptide_df))
    peptide_df = combined_peptide_df.drop_duplicates(subset=['Peptide:Protein']) # Goes from 3963 to 1440 peptides for wt v mutant dataset
    # print(len(peptide_df))
    # DF that contains all unique non-tryptic peptides across all samples
    master_peptide_df = pd.DataFrame(peptide_df, columns=['Peptide:Protein', 'Tryptic State', 'Protein ID', 'Prev AA', 'Next AA']) 

    # DF that contains all non-tryptic peptides across all samples (3963 non-tryptic peptides)
    master_peptide_df_all = pd.DataFrame(combined_peptide_df, columns=['Peptide:Protein', 'Tryptic State', 'Protein ID', 'Prev AA', 'Next AA']) 

    split_data = master_peptide_df['Peptide:Protein'].str.split(':', n=1, expand=True)
    master_peptide_df['Peptide'] = split_data[0]
    master_peptide_df['Protein'] = split_data[1]

    split_data_all = master_peptide_df_all['Peptide:Protein'].str.split(':', n=1, expand=True)
    master_peptide_df_all['Peptide'] = split_data_all[0]
    master_peptide_df_all['Protein'] = split_data_all[1]

    for i, df in enumerate(processed_dict.values()):
        sample_name = list(processed_dict.keys())[i]
        # if len(sample_name) > 2:
        master_peptide_df[sample_name] = master_peptide_df['Peptide:Protein'].isin(df['Peptide:Protein']).astype(int)
        master_peptide_df_all[sample_name] = master_peptide_df_all['Peptide:Protein'].isin(df['Peptide:Protein']).astype(int)

    # File that contains all non-tryptic peptides across samples
    master_peptide_df_all.to_csv(os.path.join(output_directory, 'all_non_tryptic_peptides.csv'), index=False)
    # Get counts of each protein associated with non-tryptic peptides (how many times a certain protein appears among the list of all non-tryptic peptides)
    # Gives us qualitative insight on potential targets of endogenous proteolysis (if a certain protein has many occurrences of non-tryptic peptides = more likely target of endogenous proteolysis)
    protein_counts = master_peptide_df_all['Protein ID'].value_counts()
    print(protein_counts)
    protein_counts.to_csv(os.path.join(output_directory, 'non-tryptic_protein_counts.csv'), index=False)


############ Protein-level output ############ 
    # Drop duplicates based on protein column to get only one row per protein, presence still based on peptide associated with this protein
    proteins_df = master_peptide_df.drop_duplicates(subset=['Protein ID'])
    proteins_df.to_csv(os.path.join(output_directory, 'non-tryptic_unique_proteins.csv'), index=False)

    return master_peptide_df



def fetch_sequence(peptide):
    url = f"https://www.uniprot.org/uniprot/{peptide}.fasta" 
    # response = requests.get(url)
    # Check if successful request and extract protein sequence 
    try:
    # if response.ok:
        response = requests.get(url, timeout=15)
        response.raise_for_status()
        lines = response.text.strip().split("\n")
        protein_sequence = "".join(lines[1:])
        return peptide, protein_sequence
    except requests.exceptions.RequestException: 
    # else:
        return peptide, None

def get_protein_sequence(peptide_df):
    peptides_id = peptide_df['Protein ID'].to_list()
    # Can alter this value to set the number of worker threads for parallel execution 
    max_workers = 50

    # Store protein sequences with respective protein IDs in dictionary
    protein_sequences = {}

    # Use ThreadPoolExecutor from concurrent.futures to launch parallel fetching of protein sequences
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        # List to store Future objects which represent the execution of the fetch_sequence method
        futures = []
        for peptide in peptides_id:
            # Submit fetch_sequence method for execution and store the Future object and then add it to the futures list
            future = executor.submit(fetch_sequence, peptide)
            futures.append(future)

        # Iterate over completed Future objects
        for future in concurrent.futures.as_completed(futures):
            # Extract the result of the fetch_sequence method from the Future object using the .result method
            peptide, sequence = future.result()
            # Add sequence to respective protein ID in dictionary
            if sequence:
                protein_sequences[peptide] = sequence

        
        # Create Protein Sequence column in df and output to file
        peptide_df['Protein Sequence'] = [protein_sequences[peptide] for peptide in peptides_id]
        # peptide_df.to_csv(f'{output_directory}/nontryp_pept_sequences.csv', index=False)
        # print(peptide_df.columns)
        return peptide_df


def peptide_seq_match(df:pd.DataFrame, output_directory):
    # match peptide with protein sequence using re 
    peptide = df['Peptide'].astype(str)
    prev_AA = df['Prev AA'].astype(str)
    tryp_state = df['Tryptic State'].astype(str)
    sequence = df['Protein Sequence'].astype(str)

    n_terminal_dict = {}
    c_terminal_dict = {}

    for peptide, prev_AA, tryp_state, sequence in zip(peptide, prev_AA, tryp_state, sequence):
        matches = list(re.finditer(peptide, sequence))
        if matches:
            match = matches[0]
            start_index = match.start()
            end_index = match.end()
    
            # Calculate indices for extracting the 5 amino acids before and after the match
            # 2nd argument acts as a boundary, so there's no values that go outside of the sequence length that can't be indexed (can't go below 0, can't go above length of sequence)
            before_start_index = max(start_index - 4, 0)
            after_end_index = min(end_index + 4, len(sequence))

            # Extract 5 amino acids before and after the match
            before_match = sequence[before_start_index:start_index]
            after_match = sequence[end_index:after_end_index]

            n_terminal = before_match + peptide[:4]
            c_terminal = peptide[-4:] + after_match
            # print(f"N-Terminus: {n_terminal}, C-Terminus: {c_terminal}")
            
            if tryp_state == 'Non-Tryptic':
                n_terminal_dict[peptide] = n_terminal
                c_terminal_dict[peptide] = c_terminal
            else: # If its semi-tryptic
                if prev_AA in ['K', 'R']:
                    c_terminal_dict[peptide] = c_terminal
                    n_terminal_dict[peptide] = pd.NA
                elif peptide.endswith('K') or peptide.endswith('R'):
                    n_terminal_dict[peptide] = n_terminal
                    c_terminal_dict[peptide] = pd.NA

    # # Add dashes to the end of sequences that are less than 10 aa long for N-terminal cleavages
    for peptide in n_terminal_dict:
        string = n_terminal_dict[peptide]
        if pd.isna(string):
            continue
        while len(string) < 8:
            if len(string) == 0:
                break
            else:
                string = "-" + string
        n_terminal_dict[peptide] = string
    # # Add dashes to beginning of sequences that are less than 10 aa long for C-terminal cleavages
    for peptide in c_terminal_dict:
        string = c_terminal_dict[peptide]
        if pd.isna(string):
            continue
        while len(string) < 8:
            if len(string) == 0:
                break
            else:
                string += "-"
        c_terminal_dict[peptide] = string
                
    df['N-terminal'] = df['Peptide'].map(n_terminal_dict)
    df['C-terminal'] = df['Peptide'].map(c_terminal_dict)

    important_cols = ['Peptide:Protein','Tryptic State','Protein ID','Prev AA','Next AA','Peptide','Protein','Protein Sequence','N-terminal','C-terminal']
    new_cols = important_cols + (df.columns.drop(important_cols).tolist())
    df = df[new_cols]
    df.to_csv(f'{output_directory}/non-tryptic_pept_sequences.csv', index=False)
    return df

    # nontryp_termini_only = nontryp_termini.dropna()

def create_termini_list(df:pd.DataFrame, output_directory):
    df_subset = df.iloc[:, 10:]
    presence_cols = list(df_subset.columns.values)
    df_N = df[['N-terminal', 'Protein', 'Protein ID', 'Tryptic State'] + presence_cols].rename(columns={'N-terminal': 'Non-Tryptic Termini'})
    df_C = df[['C-terminal', 'Protein', 'Protein ID', 'Tryptic State'] + presence_cols].rename(columns={'C-terminal': 'Non-Tryptic Termini'})

# Concatenate these two dataframes
    nontryp_termini = pd.concat([df_N, df_C])
    # nontryp_termini = pd.concat([df_subset[['N-terminal'] + presence_cols].rename(columns={'N-terminal': 'Non-Tryptic Termini'}),
    #                              df_subset[['C-terminal'] + presence_cols].rename(columns={'C-terminal': 'Non-Tryptic Termini'})])
    # nontryp_termini = pd.concat([df[['N-terminal', 'Protein ID', 'Protein']].rename(columns={'N-terminal': 'Non-Tryptic Termini'}),
    #                              df[['C-terminal', 'Protein ID', 'Protein']].rename(columns={'C-terminal': 'Non-Tryptic Termini'})])
    nontryp_termini = nontryp_termini.dropna()
    nontryp_termini.to_csv(f'{output_directory}/non-tryptic_termini.csv', index=False)
    return nontryp_termini




