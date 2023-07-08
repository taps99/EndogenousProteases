import tkinter as tk
from tkinter import ttk
from tkinter import filedialog
from tkinter import messagebox
import pandas as pd
import os
from new_functions import process_data, peptide_seq_match, get_protein_sequence, output_files
import threading
import time


output_directory = ""
# Function to open file explorer to select files to be processed
def browse_files():
    count = 0
    processed_filepaths = []

    file_paths = filedialog.askopenfilenames(filetypes=[("TSV files", "*.tsv"), ("All files", "*.*")], multiple=True)
    for file_path in file_paths:
        file_listbox.insert(tk.END, file_path)

# Function to remove selected files from list
def remove_files():
    selected_files = file_listbox.curselection()
    for index in selected_files[::-1]:
        file_listbox.delete(index)

# Function to select the directory that the output files should be output to
def select_output_directory():
    global output_directory
    output_directory = filedialog.askdirectory()
    if len(str(output_directory)) > 0:
        output_directory_label["text"] = f"Output directory: {output_directory}"
    os.makedirs(output_directory, exist_ok=True)
    return output_directory
    
# Function to setup interface elements (count of files, data structures to hold data,)
def setup():
    count = 0
    processed_filepaths = []
    listbox_entries = []

    for i, listbox_entry in enumerate(file_listbox.get(0, tk.END)):
        file_path = file_listbox.get(i)
        listbox_entries.append(listbox_entry)
        directory_name = os.path.basename(os.path.dirname(file_path))
        # file_name_without_ext = os.path.splitext(os.path.basename(file_path))[0]
        # output_file_name = f"{file_name_without_ext}_output_{directory_name}.tsv"
        sample_name = directory_name.split('.')[0]
        processed_filepaths.append(sample_name)
        
    return processed_filepaths, listbox_entries

# This gets executed when you click 'Process' button

def process_files_button():

    if not file_listbox.get(0, tk.END):
        messagebox.showerror("Error", "Please select input files before processing.")
        completion_label["text"] = f"Please select a file to process."
        return

    if not output_directory:
        messagebox.showerror("Error", "Please select output directory before processing.")
        return
    
    threading.Thread(target=processing_functions).start()
    # processing_functions()
    # peptide_seq_match(get_protein_sequence(output_files(processed_filepaths, process_data(listbox_entries)[0], process_data(listbox_entries)[1], output_directory)), output_directory)


    # args = ((output_files(processed_filepaths, processed_dfs, unprocessed_dfs)),)
    # threading.Thread(target=get_protein_sequence, args=args).start()
    # output_files(processed_filepaths, processed_dfs, unprocessed_dfs))

def processing_functions():

    processed_filepaths = setup()[0]
    listbox_entries = setup()[1]
    for i in range(0, 35):
        barnum.set(i)
        time.sleep(0.05)
    non_tryp_peptides = output_files(processed_filepaths, process_data(listbox_entries)[0], process_data(listbox_entries)[1], output_directory)
    for i in range(35, 50):
        barnum.set(i)
        time.sleep(0.05)
    for i in range(50, 81):
        barnum.set(i)
        time.sleep(0.05)
    peptide_seq_match(get_protein_sequence(non_tryp_peptides), output_directory)
    for i in range(80, 101):
        barnum.set(i)
        time.sleep(0.05)
    count = len(listbox_entries)
    if count == 1:
        completion_label["text"] = f"Successfully processed {count} file."
    else:
        completion_label["text"] = f"Successfully processed {count} files."

    output_filepaths = []
    for dirpath, dirnames, filenames in os.walk(output_directory):
        for filename in filenames:
            output_filepaths.append(os.path.join(dirpath, filename))
    paths_text = '\n'.join(output_filepaths)
    output_files_label["text"] = f"Output files:\n{paths_text}"


# Window and widgets
window = tk.Tk()
window.title("Non-Tryptic Peptide Extractor")
window.geometry("800x600")
barnum = tk.IntVar()
# Style configuration of widgets
style = ttk.Style()
window.tk.call("source", "Azure-ttk-theme-main/azure.tcl")
window.tk.call("set_theme", "dark")

style.configure("TButton", padding=5, font=('TkDefaultFont'))
style.configure("TListbox", padding=10)
style.configure("TProgressBar", thickness=50)


notebook = ttk.Notebook(window)
notebook.grid()
main_tab = tk.Frame(notebook, width=2400, height=1600)
main_tab.grid()
second_tab = tk.Frame(notebook, width=2400, height=1600)
second_tab.grid()
notebook.add(main_tab, text="Non-Tryptic Peptide Extraction")
notebook.add(second_tab, text="2nd tab")

# Define + arrange widgets in the interface
frame_browse = tk.Frame(main_tab)
frame_browse.grid(row=2, column=0, sticky=tk.W)
browse_button = ttk.Button(frame_browse, text="Browse Files", command=browse_files)
browse_button.grid(row=2, column=0, pady=10, padx=10, sticky=tk.W)

frame_remove = tk.Frame(main_tab)
frame_remove.grid(row=2, column=1, sticky=tk.W)
remove_button = ttk.Button(frame_browse, text="Remove Selected", command=remove_files)
remove_button.grid(row=2, column=1, pady=10, padx=10, sticky=tk.W)

frame_process = tk.Frame(main_tab)
# frame_process.grid(row=2, column=2, sticky=tk.W)
# frame_process.place(x=10,y=342)
process_button = ttk.Button(main_tab, text="Process", command=process_files_button, style='Accent.TButton')
process_button.grid(row=4, column=0, pady=10, padx=10, sticky=tk.W)
# process_button.place(x=10,y=342)

listbox_label = tk.Label(main_tab, text="List of selected PSM files:")
listbox_label.grid(row=0, column=0, pady=10, padx=10, sticky=tk.W)

file_listbox = tk.Listbox(main_tab, selectmode=tk.MULTIPLE, width=100, height=10)
file_listbox.grid(row=1, column=0, columnspan=2, padx=10, sticky=tk.W)

completion_label = tk.Label(main_tab, text="", fg='#008BFF', font=("TkDefaultFont", 11, "bold"), justify= tk.LEFT)
completion_label.grid(row=9, column=0, pady=10, padx=5, sticky=tk.W)

# frame_output_files = tk.Frame(window)
# frame_output_files.grid(row=6, column=0, sticky=tk.S)
# output_listbox_label = tk.Label(frame_output_files, text="Output files:")
# output_listbox_label.grid(row=6, column=0, pady=10, padx=10, sticky=tk.S)

# frame_output_files = tk.Frame(window)
# frame_output_files.grid(row=10, column=0, sticky=tk.W)
# output_files_label = tk.Label(frame_output_files, text="", anchor="w")
# output_files_label.grid(row=10, column=0, pady=10, padx=10, sticky=tk.W+tk.E)

frame_output_files = tk.Frame(main_tab)
frame_output_files.grid(row=12, column=0, sticky=tk.W)
output_files_label = tk.Label(main_tab, text="", justify= tk.LEFT)
output_files_label.grid(row=12, column=0,pady=10, padx=10, sticky='w')

frame_output_directory = tk.Frame(main_tab)
frame_output_directory.grid(row=3, column=0, sticky=tk.W)
output_directory_button = ttk.Button(frame_output_directory, text="Select Output Directory", command=select_output_directory)
output_directory_button.grid(row=3, column=0, pady=10, padx=10, sticky=tk.W)

frame_output_directory_label = tk.Frame(main_tab)
frame_output_directory_label.grid(row=3, column=1, sticky=tk.W)
output_directory_label = tk.Label(frame_output_directory, text="No output directory selected yet.")
output_directory_label.grid(row=3, column=1, pady=10, padx=10, sticky=tk.W)


loading_bar = ttk.Progressbar(main_tab, orient='horizontal', mode='determinate', length=200, variable=barnum, maximum=100.0)
loading_bar.grid(row=8, column=0, pady=10, padx=10, sticky=tk.W)

# output_listbox = tk.Listbox(.window, selectmode=tk.MULTIPLE, width=100, height=10)
# output_listbox.grid(row=6, column=0, columnspan=2, padx=10, sticky=tk.W)


window.mainloop()