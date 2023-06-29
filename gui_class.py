import tkinter as tk
from tkinter import ttk
from tkinter import *
from tkinter.ttk import *
from tkinter import filedialog
from tkinter import messagebox
import pandas as pd
import os

class GUI:
    def __init__(self, window):
        self.window = window
        self.processed_dfs = []
        self.processed_filepaths = []
        self.output_directory = None
        self.file_listbox = None
        self.completion_label = None
        self.output_listbox = None

        # Setup widgets
        self.setup_widgets()


    def setup_widgets(self):
        # window = tk.Tk()
        self.window.title("Non-Tryptic Peptide Extractor")
        self.window.geometry("800x600")

        # Style configuration
        style = ttk.Style()
        style.configure("TButton", padding=5, background="#ccc", font=('TkDefaultFont'))
        style.configure("TListbox", padding=10, background="#eee")
        style.configure('bold.TButton', padding=5, background="blue", font=('TkDefaultFont', 9,"bold"))

        # Widgets
        frame1 = tk.Frame(self.window)
        frame1.grid(row = 2, column=0, sticky=tk.W)

        frame2 = tk.Frame(self.window)
        frame2.grid(row = 3, column=0, sticky=tk.W)

        browse_button = ttk.Button(frame1, text="Browse Files", command=self.browse_files)
        browse_button.grid(row=2, column=0, pady=10,padx=10, sticky=tk.W)
        # browse_button.pack(in_=top, side=LEFT)

        remove_button = ttk.Button(frame1, text="Remove Selected", command=self.remove_files)
        remove_button.grid(row=2, column=1, pady=10, padx=10, sticky=tk.W)
        # remove_button.pack(in_=top, side=LEFT)

        process_button = ttk.Button(self.window, text="Process", command=self.process_files, style='bold.TButton')
        process_button.grid(row=2, column=1, pady=10, padx=10, sticky=tk.E)

        listbox_label = ttk.Label(self.window, text="List of selected PSM files:")
        listbox_label.grid(row=0, column=0, pady=10, padx=10, sticky=tk.W)
        # listbox_label.pack(in_=top, side=LEFT)

        self.file_listbox = tk.Listbox(self.window, selectmode=tk.MULTIPLE, width=100, height=10)
        self.file_listbox.grid(row=1, column=0, columnspan=2, padx=10, sticky=tk.W)
        # file_listbox.pack(in_=top, side=TOP, expand=TRUE)

        self.completion_label = Label(self.window, text=" ")
        self.completion_label.grid(row=7, column=0, pady=10, padx=10, sticky=tk.W)

        output_listbox_label = Label(self.window, text="Output files:")
        output_listbox_label.grid(row=5, column=0, pady=10, padx=10, sticky=tk.W)

        self.output_listbox = tk.Listbox(self.window, selectmode=tk.MULTIPLE, width=100, height=10)
        self.output_listbox.grid(row=6, column=0, columnspan=2, padx=10, sticky=tk.W)

        output_directory_button = ttk.Button(frame2, text="Select Output Directory", command=self.select_output_directory)
        output_directory_button.grid(row=3, column=0, pady=10, padx=10, sticky=tk.W)

        self.output_directory_label = Label(frame2, text="No output directory selected yet.")
        self.output_directory_label.grid(row=3, column=1, pady=10, padx=10, sticky=tk.W)


    def browse_files(self):
        file_paths = filedialog.askopenfilenames(filetypes=[("TSV files", "*.tsv"), ("All files", "*.*")], multiple=True)
        for file_path in file_paths:
            self.file_listbox.insert(tk.END, file_path)

    def remove_files(self):
        selected_files = self.file_listbox.curselection()
        for index in selected_files[::-1]:
            self.file_listbox.delete(index)

    def select_output_directory(self):
        global output_directory
        output_directory = filedialog.askdirectory()
        self.output_directory_label["text"] = "Output directory: " + output_directory
        

    def setup(self):
        global processed_dfs
        global processed_filepaths
        processed_dfs = []
        processed_filepaths = []
        count = 0

        os.makedirs(output_directory, exist_ok=True)

        for i, listbox_entry in enumerate(self.file_listbox.get(0, END)):
            file_path = self.file_listbox.get(i)
            directory_name = os.path.basename(os.path.dirname(file_path))
            # file_name_without_ext = os.path.splitext(os.path.basename(file_path))[0]
            # output_file_name = f"{file_name_without_ext}_output_{directory_name}.tsv"
            sample_name = directory_name.split('.')[0]
            processed_filepaths.append(sample_name)
            
            count += 1

        if count == 0:
            self.completion_label["text"] = f"Please select a file to process."
        elif count == 1:
            self.completion_label["text"] = f"Successfully processed {count} file."
        else:
            self.completion_label["text"] = f"Successfully processed {count} files."
        return processed_dfs, processed_filepaths


    def process_data(self, processed_filepaths):
        unprocessed_dfs = []
        count = 0
        for i, listbox_entry in enumerate(self.file_listbox.get(0, END)):
            # Read the TSV file
            count += 1
            df = pd.read_csv(listbox_entry, delimiter='\t', low_memory=False)
            unprocessed_dfs.append(df)

            # Data transformations
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
            
            processed_dfs.append(merged_df_subset)
        return processed_dfs, unprocessed_dfs


    def output_files(self, processed_filepaths, processed_dfs, unprocessed_dfs):
        # Your code to output the files here...
        processed_dict = dict(zip(processed_filepaths, processed_dfs))
        unprocessed_dict = dict(zip(processed_filepaths, unprocessed_dfs))
        # print(processed_dict_peptide)

        # This is to get the file that contains ALL peptides from the input files (includes tryptic and non-tryptic)
        combined_df = pd.concat(unprocessed_dfs, ignore_index=TRUE)
        combined_df["Peptide:Protein"] = combined_df[['Peptide', 'Protein Description']].agg(':'.join, axis=1)
        # combined_df.to_csv('./all_peptides.csv', index=False)
        combined_df.to_csv(os.path.join(output_directory, 'all_peptides.csv'), index=False)


        # This is to get the unique proteins from the entire dataset (includes tryptic and non-tryptic)
        combined_df_subset = combined_df[['Peptide:Protein', 'Peptide', 'Prev AA', 'Next AA', 'Protein Description', 'Protein ID']].drop_duplicates(subset='Protein ID')
        for i, df in enumerate(unprocessed_dict.values()):
            sample_name = list(unprocessed_dict.keys())[i]
            for protein in combined_df_subset:
                combined_df_subset[sample_name] = combined_df_subset['Protein ID'].isin(df['Protein ID']).astype(int)

        # combined_df_subset.to_csv('all_unique_proteins.csv', index=False)
        combined_df_subset.to_csv(os.path.join(output_directory, 'all_unique_proteins.csv'), index=False)
        # result_df = pd.DataFrame(columns=['Peptide:Protein'])


        # This part involves the actual output of non-tryptic peptides + relevant information
        unique_peptide_dfs=[]
        for df in processed_dfs: # GETS UNIQUE PEPTIDES FROM EACH SAMPLE
            unique_peptides = df.drop_duplicates(subset=['Peptide:Protein'])
            unique_peptide_dfs.append(unique_peptides) 
            # print(unique_peptides)


        combined_peptide_df = pd.concat(unique_peptide_dfs) # CONCAT THE UNIQUE PEPTIDES FROM EACH SAMPLE INTO ONE LARGE DF
        # master_df = combined_df_final['Peptide:Protein'].unique() # THEN ONLY KEEP UNIQUE ONES, WHICH MEANS THIS ONE CONTAINS ALL UNIQUE PEPTIDES FROM ALL SAMPLES
        peptide_df = combined_peptide_df.drop_duplicates(subset=['Peptide:Protein'])
        master_peptide_df = pd.DataFrame(peptide_df, columns=['Peptide:Protein', 'Tryptic State', 'Protein ID', 'Prev AA', 'Next AA']) # DF that contains all unique peptides across all samples

    
        split_data = master_peptide_df['Peptide:Protein'].str.split(':', n=1, expand=True)
        master_peptide_df['Peptide'] = split_data[0]
        master_peptide_df['Protein'] = split_data[1]


        for i, df in enumerate(processed_dict.values()):
            sample_name = list(processed_dict.keys())[i]
            for peptide in master_peptide_df:
                master_peptide_df[sample_name] = master_peptide_df['Peptide:Protein'].isin(df['Peptide:Protein']).astype(int)
        
        # File that contains all non-tryptic peptides across samples
        # master_peptide_df.to_csv('./non_tryptic_peptides.csv', index=False)
        master_peptide_df.to_csv(os.path.join(output_directory, 'non_tryptic_peptides.csv'), index=False)



    ############ Protein-level output ############ 
        # Drop duplicates based on protein column to get only one row per protein, presence still based on peptide associated with this protein
        proteins_df = master_peptide_df.drop_duplicates(subset=['Protein ID'])
        # proteins_df.to_csv('./non_tryptic_proteins.csv', index=False)
        proteins_df.to_csv(os.path.join(output_directory, 'non_tryptic_proteins.csv'), index=False)


        output_filepaths = ['./all_peptides.csv', './all_unique_proteins.csv', './non_tryptic_peptides.csv', './protein_output.csv']
        for filepath in output_filepaths:
            df = pd.read_csv(filepath, low_memory=False)
            peptide_count = df['Peptide'].nunique() # Modify this line if your column name differs
            self.output_listbox.insert(tk.END, f'File: {filepath}, Number of peptides: {peptide_count}')

        return processed_dfs


    def process_files(self):
        self.processed_dfs, self.processed_filepaths = self.setup()
        self.processed_dfs, unprocessed_dfs = self.process_data(self.processed_filepaths)
        self.output_files(processed_filepaths, processed_dfs, unprocessed_dfs)



def main():
    # Initialize window
    window = tk.Tk()
    gui = GUI(window)
    window.mainloop()


if __name__ == "__main__":
    main()
