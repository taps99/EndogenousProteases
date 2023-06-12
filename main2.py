import tkinter as tk
from tkinter import ttk
from tkinter import *
from tkinter.ttk import *
from tkinter import filedialog
import pandas as pd
import os
import numpy as np


def browse_files():
    file_paths = filedialog.askopenfilenames(filetypes=[("TSV files", "*.tsv"), ("All files", "*.*")])
    for file_path in file_paths:
        file_listbox.insert(tk.END, file_path)

def remove_files():
    selected_files = file_listbox.curselection()
    for index in selected_files[::-1]:
        file_listbox.delete(index)

def process_files():
    count = 0
    global processed_dfs
    global processed_filepaths
    global sample_name
    processed_dfs = []
    processed_filepaths = []
    # processed_dict = dict()
    

    # processed_file_listbox.delete(0, tk.END)
    for i, listbox_entry in enumerate(file_listbox.get(0, END)):
        file_path = file_listbox.get(i)
        directory_name = os.path.basename(os.path.dirname(file_path))
        file_name_without_ext = os.path.splitext(os.path.basename(file_path))[0]

        # Save the processed DataFrame to an output file
        output_file_name = f"{file_name_without_ext}_output_{directory_name}.tsv"
        sample_name = directory_name.split('.')[0]
        
        processed_filepaths.append(sample_name)
        count +=1
    # for file in processed_filepaths:
    #     # processed_file_listbox.insert(0, file) # Shows list of processed files in 2nd listbox
    #     sample_name = file.split('.')[0]
    #     count +=1

    # for i in range(file_listbox.size()):
    for i, listbox_entry in enumerate(file_listbox.get(0, END)):
        # Read the TSV file
        df = pd.read_csv(listbox_entry, delimiter='\t', low_memory=False)
        
        # Filter "Prev AA" column to exclude peptides that have K or R preceding them (This results in semi-tryptic peptides)
        mask_prev = (df["Prev AA"] != "K") & (df["Prev AA"] != "R")
        df_prevnoKR = df[mask_prev]

        # From this, we can filter out peptides that terminate with K or R as well
        mask_end = ~(df["Peptide"].str.endswith("K") | df["Peptide"].str.endswith("R"))
        df_endnoKR = df[mask_end]

        # Concatenate the df containing semi-tryptic peptides without preceding K/R's, and the df containing semi-tryptic peptides that don't terminate with K/R's
        merged_df = pd.concat([df_prevnoKR, df_endnoKR]).drop_duplicates().reset_index(drop=True)
        

        prev_aa_mask = merged_df['Prev AA'].str.contains('K|R')
        peptide_end_mask = merged_df['Peptide'].str.endswith('K') | merged_df['Peptide'].str.endswith('R')

        merged_df.loc[prev_aa_mask | peptide_end_mask, 'Tryptic State'] = 'Semi-Tryptic'
        merged_df.loc[~(prev_aa_mask | peptide_end_mask), 'Tryptic State'] = 'Non-Tryptic'

        # Subset columns of interest that are useful to import into Perseus for further analysis
        # merged_df_subset = merged_df[["Spectrum", "Peptide", "Prev AA", "Intensity", "Protein ID", "Entry Name", "Gene", "Mapped Proteins", "Protein Description", "Tryptic State"]]
        merged_df["Peptide:Protein"] = merged_df[['Peptide', 'Protein Description']].agg(':'.join, axis=1)
        merged_df_subset = merged_df[['Peptide:Protein', 'Tryptic State', 'Peptide', 'Protein Description', 'Protein ID']]
        
        processed_dfs.append(merged_df_subset)
        



######## Peptide-level output #######
    processed_dict = dict(zip(processed_filepaths, processed_dfs))
    # print(processed_dict_peptide)

    combined_df = pd.concat(processed_dfs, ignore_index=TRUE)
    result_df = pd.DataFrame(columns=['Peptide:Protein'])


    unique_peptide_dfs=[]
    for df in processed_dfs: # GETS UNIQUE PEPTIDES FROM EACH SAMPLE
        unique_peptides = df.drop_duplicates(subset=['Peptide:Protein'])
        unique_peptide_dfs.append(unique_peptides) 
        # print(unique_peptides)


    combined_peptide_df = pd.concat(unique_peptide_dfs) # CONCAT THE UNIQUE PEPTIDES FROM EACH SAMPLE INTO ONE LARGE DF
    # master_df = combined_df_final['Peptide:Protein'].unique() # THEN ONLY KEEP UNIQUE ONES, WHICH MEANS THIS ONE CONTAINS ALL UNIQUE PEPTIDES FROM ALL SAMPLES
    peptide_df = combined_peptide_df.drop_duplicates(subset=['Peptide:Protein'])
    master_peptide_df = pd.DataFrame(peptide_df, columns=['Peptide:Protein', 'Tryptic State', 'Protein ID']) # DF that contains all unique peptides across all samples
    
    print(len(master_peptide_df))
 
    split_data = master_peptide_df['Peptide:Protein'].str.split(':', n=1, expand=True)
    master_peptide_df['Peptide'] = split_data[0]
    master_peptide_df['Protein'] = split_data[1]

    # master_df['Tryptic State'] = master_df['Peptide:Protein'].map(combined_df.set_index('Peptide:Protein')['Tryptic State'])

    for i, df in enumerate(processed_dict.values()):
        sample_name = list(processed_dict.keys())[i]
        for peptide in master_peptide_df:
            master_peptide_df[sample_name] = master_peptide_df['Peptide:Protein'].isin(df['Peptide:Protein']).astype(int)
    
    master_peptide_df.to_csv('./peptide_output.csv', index=False)



############ Protein-level output ############ 

    unique_protein_dfs = []
    for df in processed_dfs:
        unique_proteins = df.drop_duplicates(subset=['Protein Description'])
        unique_protein_dfs.append(unique_proteins) # Generates list of dataframes containing only the unique proteins from each sample
    
    combined_protein_df = pd.concat(unique_protein_dfs) # CONCAT THE UNIQUE PROTEINS FROM EACH SAMPLE INTO ONE LARGE DF
    protein_df = combined_protein_df.drop_duplicates(subset=['Protein Description'])
    master_protein_df = pd.DataFrame(protein_df, columns=['Protein Description', 'Tryptic State', 'Protein ID']) # DF that contains all unique peptides across all samples

    print(len(master_protein_df))

    for i, df in enumerate(processed_dict.values()):
        sample_name = list(processed_dict.keys())[i]
        for protein in master_protein_df:
            master_protein_df[sample_name] = master_protein_df['Protein Description'].isin(df['Protein Description']).astype(int)
    
    master_protein_df.to_csv('./protein_output.csv', index=False)


    if count == 0:
        completion_label["text"] = f"Please select a file to process."
    elif count == 1:
        completion_label["text"] = f"Successfully processed {count} file."
    else:
        completion_label["text"] = f"Successfully processed {count} files."


    # merged_df_subset.to_csv('./{output_file_name}', sep='\t', index=False)

    return processed_dfs  


# def merge_files():
#     unique_list = []
#     bleh = []
#     combined_df = pd.concat(processed_dfs, ignore_index=TRUE)

#     result_df = pd.DataFrame(columns=['Peptide:Protein'])


#     unique_list=[]
#     for df in processed_dfs: # GETS UNIQUE PEPTIDES FROM EACH SAMPLE
#         # df_all = df1.merge(df.drop_duplicates(), on=['Peptide:Protein', 'Tryptic State'], how='left', indicator=True)
#         unique_peptides = df['Peptide:Protein'].unique() # Generates a numpy array of all unique peptides found in a given sample
#         unique_list.append(unique_peptides) # List of numpy arrays, which should contain all unique peptides for all samples

#     bleh = []
#     for arr in unique_list: # CONVERT ARRAYS TO DFs
#         unique_df = pd.DataFrame(arr, columns=['Peptide:Protein'])
#         bleh.append(unique_df)
#     # print(len(bleh))

#     combined_df_final = pd.concat(bleh) # CONCAT THE UNIQUE PEPTIDES FROM EACH SAMPLE INTO ONE LARGE DF
#     master_df = combined_df_final['Peptide:Protein'].unique() # THEN ONLY KEEP UNIQUE ONES, WHICH MEANS THIS ONE CONTAINS ALL UNIQUE PEPTIDES FROM ALL SAMPLES
#     master_df = pd.DataFrame(master_df, columns=['Peptide:Protein']) # DF that contains all unique peptides across all samples
    
#     print(len(master_df))
 

#     for i, df in enumerate(processed_dfs, start=1):
#         sample_name = f'Sample {i}'

#         for peptide in master_df:
#             master_df[sample_name] = master_df['Peptide:Protein'].isin(df['Peptide:Protein']).astype(int)
    
#     # print(master_df.columns[1:])
    

#     master_df.to_csv('./final_output.csv', index=False)
        





### Customize window and buttons 
# Initialize window
window = tk.Tk()
window.title("Non-Tryptic Peptide Extractor")
window.geometry("700x400")

# Style configuration
style = ttk.Style()
# style.configure("TLabel", font=(None, 12))
style.configure("TButton", padding=6, background="#ccc")
style.configure("TListbox", padding=10, background="#eee")

# top = ttk.Frame(window)
# top.pack()
# bottom = ttk.Frame(window)
# bottom.pack(side=BOTTOM)

# Widgets
browse_button = ttk.Button(window, text="Browse Files", command=browse_files)
browse_button.grid(row=2, column=0, pady=10, padx=10, sticky=tk.W)
# browse_button.pack(in_=top, side=LEFT)

remove_button = ttk.Button(window, text="Remove Selected", command=remove_files)
remove_button.grid(row=2, column=0, pady=10, padx=10, sticky=tk.S)
# remove_button.pack(in_=top, side=LEFT)

process_button = ttk.Button(window, text="Process", command=process_files)
process_button.grid(row=2, column=1, pady=10, sticky=tk.E)

listbox_label = ttk.Label(window, text="List of selected PSM files:")
listbox_label.grid(row=0, column=0, pady=10, padx=10, sticky=tk.W)
# listbox_label.pack(in_=top, side=LEFT)

file_listbox = tk.Listbox(window, selectmode=tk.MULTIPLE, width=100, height=10)
file_listbox.grid(row=1, column=0, columnspan=2, padx=10, sticky=tk.W)
# file_listbox.pack(in_=top, side=TOP, expand=TRUE)

completion_label = Label(window, text=" ")
completion_label.grid(row=3, column=0, pady=10, padx=10, sticky=tk.W)

# pfl_label = Label(window, text="List of processed files:")
# pfl_label.grid(row=4, column=0, pady=10, padx=10, sticky=tk.W)

# processed_file_listbox = tk.Listbox(window, selectmode=tk.MULTIPLE, width=100, height=10)
# processed_file_listbox.grid(row=5, column=0, columnspan=2, padx=10, sticky=tk.W)

# merge_button = ttk.Button(window, text="Merge", command=merge_files)
# merge_button.grid(row=6, column=1, pady=10, sticky=tk.E)

# process_button.pack(in_=top, side=LEFT)

# count = 0
# text_f = f"Processed {count} files."
# output_text = ttk.Label(window, text=text_f)
# output_text.pack(pady=10)

window.mainloop()