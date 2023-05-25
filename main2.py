import tkinter as tk
from tkinter import ttk
from tkinter import *
from tkinter.ttk import *
from tkinter import filedialog
import pandas as pd

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
    processed_files = []

    processed_file_listbox.delete(0, tk.END)
    for i in range (file_listbox.size()):
        file_path = file_listbox.get(i)

        # Read the TSV file
        df = pd.read_csv(file_path, delimiter='\t', low_memory=False)
        
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
        merged_df_subset = merged_df[["Spectrum", "Peptide", "Prev AA", "Intensity", "Protein ID", "Entry Name", "Gene", "Mapped Proteins", "Protein Description", "Tryptic State"]]

        # Adding a column to show whether a peptide is semi-tryptic or fully non-tryptic
        # merged_df_subset['Tryptic State'] = merged_df_subset.apply(lambda row: 'Semi-Tryptic' if row['Prev AA'] == 'K' or row['Prev AA'] == 'R' or row['Peptide'].endswith("K") or row['Peptide'].endswith("R") else 'Non-Tryptic', axis=1)

        # semi_mask = ((merged_df_subset["Prev AA"] == "K" | merged_df_subset["Prev AA"] == "R" ) | (merged_df_subset["Peptide"].str.endswith("K") | merged_df_subset["Peptide"].str.endswith("R")))
        # nontryp_mask = ((merged_df_subset["Prev AA"] != "K" | merged_df_subset["Prev AA"] != "R" ) & ~(merged_df_subset["Peptide"].str.endswith("K") | merged_df_subset["Peptide"].str.endswith("R")))

        # merged_df_subset.loc[semi_mask, "Tryptic State"] = "Semi-Tryptic"
        # merged_df_subset.loc[nontryp_mask, "Tryptic State"] = "Non-Tryptic"


        # Save the processed DataFrame to an output file
        output_file_path = file_path.replace('.tsv', '_output.tsv')
        # print(output_file_path)
        merged_df_subset.to_csv(output_file_path, sep='\t', index=False)

        processed_files.append(output_file_path)
        
    for file in processed_files:
        processed_file_listbox.insert(0, file)
        

        count +=1
    if count == 0:
        completion_label["text"] = f"Please select a file to process."
    elif count == 1:
        completion_label["text"] = f"Successfully processed {count} file."
    else:
        completion_label["text"] = f"Successfully processed {count} files."

# def merge_files():


### Customize window and buttons 
# Initialize window
window = tk.Tk()
window.title("Non-Tryptic Peptide Extractor")
window.geometry("900x700")

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

pfl_label = Label(window, text="List of processed files:")
pfl_label.grid(row=4, column=0, pady=10, padx=10, sticky=tk.W)

processed_file_listbox = tk.Listbox(window, selectmode=tk.MULTIPLE, width=100, height=10)
processed_file_listbox.grid(row=5, column=0, columnspan=2, padx=10, sticky=tk.W)

merge_button = ttk.Button(window, text="Merge", command=process_files)
merge_button.grid(row=6, column=1, pady=10, sticky=tk.E)

# process_button.pack(in_=top, side=LEFT)

# count = 0
# text_f = f"Processed {count} files."
# output_text = ttk.Label(window, text=text_f)
# output_text.pack(pady=10)

window.mainloop()