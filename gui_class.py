import tkinter as tk
from tkinter import ttk
from tkinter import filedialog
from tkinter import messagebox
import pandas as pd
import os
import requests
import concurrent.futures
import re
from tenacity import retry, stop_after_attempt, wait_fixed
import threading




# Create a class for the GUI instance
class GUI:
    # Initialize variables + window for interface
    def __init__(self, window):
        self.window = window
        self.processed_dfs = []
        self.processed_filepaths = []
        self.output_directory = None
        self.file_listbox = None
        self.completion_label = None
        # self.output_listbox = None
        self.output_files_label = None
        self.loading_bar = None
        self.listbox_entries = []

        # Setup widgets
        self.setup_widgets()

    def setup_widgets(self):
        self.window.title("Non-Tryptic Peptide Extractor")
        self.window.geometry("800x600")
        
        # Style configuration of widgets
        style = ttk.Style()
        self.window.tk.call("source", "Azure-ttk-theme-main/azure.tcl")
        self.window.tk.call("set_theme", "dark")
       
        style.configure("TButton", padding=5, font=('TkDefaultFont'))
        style.configure("TListbox", padding=10)
        style.configure('bold.TButton', padding=5, font=('TkDefaultFont', 9,"bold"), background = 'lightblue')
        # style.configure('Heading.Label',)

        # Define + arrange widgets in the interface
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

        process_button = ttk.Button(self.window, text="Process", command=self.process_files, style='Accent.TButton')
        process_button.grid(row=2, column=1, pady=10, padx=10, sticky=tk.E)

        listbox_label = tk.Label(self.window, text="List of selected PSM files:")
        listbox_label.grid(row=0, column=0, pady=10, padx=10, sticky=tk.W)
        # listbox_label.pack(in_=top, side=LEFT)

        self.file_listbox = tk.Listbox(self.window, selectmode=tk.MULTIPLE, width=100, height=10)
        self.file_listbox.grid(row=1, column=0, columnspan=2, padx=10, sticky=tk.W)
        # file_listbox.pack(in_=top, side=TOP, expand=TRUE)

        self.completion_label = tk.Label(self.window, text="")
        self.completion_label.grid(row=7, column=0, pady=10, padx=10, sticky=tk.W)

        self.output_files_label = tk.Label(self.window, text="")
        self.output_files_label.grid(row=6, column=0, pady=10, padx=10, sticky=tk.W)

        output_listbox_label = tk.Label(self.window, text="Output files:")
        output_listbox_label.grid(row=5, column=0, pady=10, padx=10, sticky=tk.W)

        # self.output_listbox = tk.Listbox(self.window, selectmode=tk.MULTIPLE, width=100, height=10)
        # self.output_listbox.grid(row=6, column=0, columnspan=2, padx=10, sticky=tk.W)

        output_directory_button = ttk.Button(frame2, text="Select Output Directory", command=self.select_output_directory)
        output_directory_button.grid(row=3, column=0, pady=10, padx=10, sticky=tk.W)

        self.output_directory_label = tk.Label(frame2, text="No output directory selected yet.")
        self.output_directory_label.grid(row=3, column=1, pady=10, padx=10, sticky=tk.W)

        self.loading_bar = ttk.Progressbar(self.window, orient='horizontal', mode='determinate', length=100)
        self.loading_bar.grid(row=8, column=0, pady=10, padx=10, sticky=tk.W)

    # Function to open file explorer to select files to be processed
    def browse_files(self):
        file_paths = filedialog.askopenfilenames(filetypes=[("TSV files", "*.tsv"), ("All files", "*.*")], multiple=True)
        for file_path in file_paths:
            self.file_listbox.insert(tk.END, file_path)

    # Function to remove selected files from list
    def remove_files(self):
        selected_files = self.file_listbox.curselection()
        for index in selected_files[::-1]:
            self.file_listbox.delete(index)

    # Function to select the directory that the output files should be output to
    def select_output_directory(self):
        # global output_directory
        self.output_directory = filedialog.askdirectory()
        self.output_directory_label["text"] = f"Output directory: {self.output_directory}"
        
    # Function to setup interface elements (count of files, data structures to hold data,)
    def setup(self):
        count = 0
        os.makedirs(self.output_directory, exist_ok=True)

        for i, listbox_entry in enumerate(self.file_listbox.get(0, tk.END)):
            file_path = self.file_listbox.get(i)
            self.listbox_entries.append(listbox_entry)
            directory_name = os.path.basename(os.path.dirname(file_path))
            # file_name_without_ext = os.path.splitext(os.path.basename(file_path))[0]
            # output_file_name = f"{file_name_without_ext}_output_{directory_name}.tsv"
            sample_name = directory_name.split('.')[0]
            self.processed_filepaths.append(sample_name)
            
            count += 1

        if count == 0:
            self.completion_label["text"] = f"Please select a file to process."
        elif count == 1:
            self.completion_label["text"] = f"Successfully processed {count} file."
        else:
            self.completion_label["text"] = f"Successfully processed {count} files."
        return self.processed_filepaths



    def process_data(self, processed_filepaths):
        self.unprocessed_dfs = []
        count = 0

        # for i, listbox_entry in enumerate(self.file_listbox.get(0, tk.END)):
        for entry in self.listbox_entries:
            # Read the TSV file
            df = pd.read_csv(entry, delimiter='\t', low_memory=False)
            self.unprocessed_dfs.append(df)

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
            
            self.processed_dfs.append(merged_df_subset)
        return self.processed_dfs, self.unprocessed_dfs


    def output_files(self, processed_filepaths, processed_dfs, unprocessed_dfs):

        processed_dict = dict(zip(processed_filepaths, processed_dfs))
        unprocessed_dict = dict(zip(processed_filepaths, unprocessed_dfs))


        # This is to get the file that contains ALL peptides from the input files (includes tryptic and non-tryptic)
        # combined_df = pd.concat(unprocessed_dfs, ignore_index=True)
        # combined_df["Peptide:Protein"] = combined_df[['Peptide', 'Protein Description']].agg(':'.join, axis=1)
        # # combined_df.to_csv(os.path.join(self.output_directory, 'all_peptides.csv'), index=False)


        # # This is to get the unique proteins from the entire dataset (includes tryptic and non-tryptic)
        # combined_df_subset = combined_df[['Peptide:Protein', 'Peptide', 'Prev AA', 'Next AA', 'Protein Description', 'Protein ID']].drop_duplicates(subset='Protein ID')
        # for i, df in enumerate(unprocessed_dict.values()):
        #     sample_name = list(unprocessed_dict.keys())[i]
        #     for protein in combined_df_subset:
        #         combined_df_subset[sample_name] = combined_df_subset['Protein ID'].isin(df['Protein ID']).astype(int)

        # combined_df_subset.to_csv(os.path.join(self.output_directory, 'all_unique_proteins.csv'), index=False)


        # This part involves the actual output of non-tryptic peptides + relevant information
        unique_peptide_dfs=[]
        for df in processed_dfs: # GETS UNIQUE PEPTIDES FROM EACH SAMPLE
            unique_peptides = df.drop_duplicates(subset=['Peptide:Protein'])
            unique_peptide_dfs.append(unique_peptides) 


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
        # master_peptide_df.to_csv(os.path.join(self.output_directory, 'non_tryptic_peptides.csv'), index=False)


    ############ Protein-level output ############ 
        # Drop duplicates based on protein column to get only one row per protein, presence still based on peptide associated with this protein
        proteins_df = master_peptide_df.drop_duplicates(subset=['Protein ID'])
        proteins_df.to_csv(os.path.join(self.output_directory, 'non_tryptic_proteins.csv'), index=False)

        return master_peptide_df

        # output_filepaths = [f'{self.output_directory}/all_peptides.csv', f'{self.output_directory}/all_unique_proteins.csv', f'{self.output_directory}/non_tryptic_peptides.csv', f'{self.output_directory}/non_tryptic_proteins.csv']
        # for filepath in output_filepaths:
        #     df = pd.read_csv(filepath, low_memory=False)
        #     peptide_count = df['Peptide'].nunique() # Modify this line if your column name differs
        #     self.output_listbox.insert(tk.END, f'File: {filepath}, Number of peptides: {peptide_count}')

    def fetch_sequence(self, peptide):
        url = f"https://www.uniprot.org/uniprot/{peptide}.fasta" 
        # response = requests.get(url)
        # Check if successful request and extract protein sequence 
        try:
        # if response.ok:
            response = requests.get(url, timeout=10)
            response.raise_for_status()
            lines = response.text.strip().split("\n")
            protein_sequence = "".join(lines[1:])
            return peptide, protein_sequence
        except requests.exceptions.RequestException: 
        # else:
            return peptide, None
    
    def get_protein_sequence(self, peptide_df):
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
                future = executor.submit(self.fetch_sequence, peptide)
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
            peptide_df.to_csv(f'{self.output_directory}/nontryp_pept_sequences.csv', index=False)
            # print(peptide_df.columns)
            return peptide_df


    def peptide_seq_match(self, df):
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
                # print(f"Match found for peptide {peptide} at index {start_index} to {end_index} in sequence {sequence}")
        
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
                        n_terminal_dict[peptide] = ""
                    elif peptide.endswith('K') or peptide.endswith('R'):
                        n_terminal_dict[peptide] = n_terminal
                        c_terminal_dict[peptide] = ""
  
        # # Add dashes to the end of sequences that are less than 10 aa long for N-terminal cleavages
        for peptide in n_terminal_dict:
            string = n_terminal_dict[peptide]
            while len(string) < 8:
                if len(string) == 0:
                    break
                else:
                    string = "-" + string
            n_terminal_dict[peptide] = string
        # # Add dashes to beginning of sequences that are less than 10 aa long for C-terminal cleavages
        for peptide in c_terminal_dict:
            string = c_terminal_dict[peptide]
            while len(string) < 8:
                if len(string) == 0:
                    break
                else:
                    string += "-"
            c_terminal_dict[peptide] = string
                    
        df['N-terminal'] = df['Peptide'].map(n_terminal_dict)
        df['C-terminal'] = df['Peptide'].map(c_terminal_dict)
        df.to_csv(f'{self.output_directory}/nontryp_pept_sequences.csv', index=False)

        return df
    

    # This gets executed when you click 'Process' button
    def process_files(self):
        self.processed_filepaths = self.setup()
        self.processed_dfs, unprocessed_dfs = self.process_data(self.processed_filepaths)
        self.peptide_seq_match(self.get_protein_sequence(self.output_files(self.processed_filepaths, self.processed_dfs, unprocessed_dfs)))
        # args = ((self.output_files(self.processed_filepaths, self.processed_dfs, unprocessed_dfs)),)
        # threading.Thread(target=self.get_protein_sequence, args=args).start()
        # self.output_files(self.processed_filepaths, self.processed_dfs, unprocessed_dfs))


def main():
    # Initialize new window, create instance of GUI class
    window = tk.Tk()
    gui = GUI(window)
    window.mainloop()


# executes main() if script is being run directly, doesn't execute if script is imported
if __name__ == "__main__":
    main()
