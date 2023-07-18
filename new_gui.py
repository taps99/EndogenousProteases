import tkinter as tk
from tkinter import ttk
from tkinter import filedialog
from tkinter import messagebox
from tkinter import simpledialog
import pandas as pd
import os
import sys
from new_functions import process_data, peptide_seq_match, get_protein_sequence, output_files, create_termini_list
import threading
import time
from matching import substrate_processing, exact_match, fuzzy_match

# def main():
output_directory = ""
organisms = []
substrate_df = None
termini_list = None

# Function to open file explorer to select files to be processed
def browse_files():
    file_paths = filedialog.askopenfilenames(filetypes=[("TSV files", "*.tsv"), ("All files", "*.*")], multiple=True)
    for file_path in file_paths:
        # file_listbox.insert(tk.END, file_path)
        treeview.insert("", tk.END, values=(file_path, ""))
    switch_button_state()

# Function to remove selected files from list
def remove_files():
    # selected_files = file_listbox.curselection()
    # for index in selected_files[::-1]:
    #     file_listbox.delete(index)
    selected_files = treeview.selection()
    for selected in selected_files:
        treeview.delete(selected)
    switch_button_state()


# Function to select the directory that the output files should be output to
def select_output_directory():
    global output_directory
    output_directory = filedialog.askdirectory()
    if len(str(output_directory)) > 0:
        output_directory_label["text"] = f"Output directory: {output_directory}"
    else: 
        output_directory_label["text"] = f"No output directory selected."
    os.makedirs(output_directory, exist_ok=True)
    switch_button_state()
    return output_directory
    
# Function to setup interface elements (count of files, data structures to hold data,)
def setup():
    count = 0
    processed_filepaths = []
    listbox_entries = []

    group_counts = {} 
    for value in treeview.get_children():
        file_path, group = treeview.item(value, "values")
        listbox_entries.append(file_path)
        # if not group:  
        #     messagebox.showerror("Error", "Please select a group for input files.")
        #     return None, None
        # Check if the group is already in group_counts
        if group in group_counts:
            # If so, increment the count and update the group name with the new count
            group_counts[group] += 1
            group = f"{group}_{group_counts[group]}"
        else:
            # If not, add it to dictionary with a count of 1 (indicating the first item)
            group_counts[group] = 1
            group = f"{group}_1"  # Append "_1" to the first item of each group
        processed_filepaths.append(group)

    return processed_filepaths, listbox_entries 

def group_files():

    if len(treeview.selection()) > 0:
        group_name = simpledialog.askstring("Input", "Enter the group name:")
        if group_name:
            selected_items = treeview.selection()
            for item in selected_items:
                file_path, _ = treeview.item(item, "values")
                treeview.item(item, values=(file_path, group_name))
    else:
        messagebox.showerror("Error", "No files were selected.")


# This gets executed when you click 'Process' button
def process_files_button():
    completion_label["text"] = ""
    output_files_label["text"] = ""

    if not treeview.get_children():
        messagebox.showerror("Error", "Please select input files before processing.")
        completion_label["text"] = f"Please select a file to process."
        return
    if not output_directory:
        messagebox.showerror("Error", "Please select output directory before processing.")
        completion_label["text"] = f"Please select an output directory."
        return
    ungrouped_files = [item for item in treeview.get_children() if not treeview.item(item, "values")[1]]
    if ungrouped_files:
        proceed = messagebox.askyesno("Warning", "Presence/Absence data will not be generated for ungrouped files. Would you like to proceed anyway?")
        if not proceed:
            return
    process_thread = threading.Thread(target=processing_functions)
    process_thread.start()

def switch_button_state():
    # if file_listbox.size() > 0 and output_directory:
    if len(treeview.get_children()) > 0 and output_directory:
        process_button['style'] = 'Accent.TButton'
    else:
        process_button['style'] = 'TButton'

def processing_functions():
    global organisms
    global substrate_df
    global termini_list
    processed_filepaths = setup()[0]
    listbox_entries = setup()[1]
    completion_label["text"] = "Initializing extraction..."
    for i in range(0, 35):
        barnum.set(i)
        time.sleep(0.05)
    try:
        processed_dfs, unprocessed_dfs = process_data(listbox_entries)
    except Exception as e:
        messagebox.showerror("Error", f"An unexpected error occurred while reading the files {listbox_entries}: {str(e)}")
        completion_label["text"] = "Process aborted. Please try again."
        barnum.set(0)
        return
    try:
        non_tryp_peptides = output_files(processed_filepaths, processed_dfs, unprocessed_dfs, output_directory)
    except Exception as e:
        messagebox.showerror("Error", f"An unexpected error occurred while generating the output files: {str(e)}")
        completion_label["text"] = "Process aborted. Please try again."
        barnum.set(0)
        return
    for i in range(35, 50):
        barnum.set(i)
        time.sleep(0.05)
    completion_label["text"] = "Processing peptides..."
    for i in range(50, 81):
        barnum.set(i)
        time.sleep(0.05)
    try:
        termini_list = create_termini_list(peptide_seq_match(get_protein_sequence(non_tryp_peptides), output_directory), output_directory)
    except Exception as e:
        messagebox.showerror("Error", f"An unexpected error occurred while acquiring sequences from UniProt: {str(e)}")
        completion_label["text"] = "Process aborted. Please try again."
        barnum.set(0)
        return
    completion_label["text"] = "Finalizing..."

    # Check if the application is running as a script or a packaged application
    if getattr(sys, 'frozen', False):
    # The application is running as a bundled executable
        application_path = sys._MEIPASS
    else:
    # The application is running as a script
        application_path = os.path.dirname(os.path.abspath(__file__))
    substrate_path = os.path.join(application_path, 'substrate.csv')
    organisms, substrate_df = substrate_processing(substrate_path)
    window.after(0, update_organism_menu, organisms)


    for i in range(80, 101):
        barnum.set(i)
        time.sleep(0.05)

    count = len(listbox_entries)
    if count == 1:
        completion_label["text"] = f"Successfully processed {count} file."
    else:
        completion_label["text"] = f"Successfully processed {count} files."

    output_filepaths = [f'{output_directory}/all_peptides.csv', f'{output_directory}/all_non_tryptic_peptides.csv', f'{output_directory}/all_unique_proteins.csv', f'{output_directory}/non-tryp_pept_sequences.csv', f'{output_directory}/non-tryptic_unique_proteins.csv', f'{output_directory}/non-tryp_termini.csv']
    # output_filepaths = [os.path.join(output_directory,f) for (dirpath, dirnames, filenames) in os.walk(output_directory) for f in filenames]
    # for dirpath, dirnames, filenames in os.walk(output_directory):
    #     for filename in filenames:
    #         output_filepaths.append(os.path.join(dirpath, filename))
    paths_text = '\n'.join(output_filepaths)
    output_files_label["text"] = f"Output files:\n{paths_text}"
    termini_count = len(termini_list)
    info_label["text"] = f"Protease matching for {termini_count} non-tryptic termini in:\n{output_filepaths[-1]}"
    match_method_menu["state"] = tk.NORMAL
    organism_menu["state"] = tk.NORMAL
    protease_button["state"] = tk.NORMAL
    protease_button["style"] = 'Accent.TButton'
    return organisms, substrate_df, termini_list


def protease_match_button():
    search_term = update_organism_menu().capitalize()
    method = update_method_menu()
    if not search_term or not method:
        messagebox.showerror("Error", "Please select a species and a matching method for protease matching.")
        return

    if not search_term:
        messagebox.showerror("Error", "Please select a species for protease matching.")
        return
    if not method:
        messagebox.showerror("Error", "Please select a matching method for protease matching.")
        return
    
    process_thread_2 = threading.Thread(target=protease_match)
    process_thread_2.start()


def update_organism_menu(*args):
    search_term = organisms_var.get().lower()
    organism_menu['values'] = [org for org in organisms if search_term in org.lower()]
    return search_term

def update_method_menu(*args):
    method = method_var.get()
    return method



tree = None
def protease_match():
    global substrate_df
    global termini_list
    global tree

    search_term = update_organism_menu().capitalize()
    method = update_method_menu()
    # print(f'Species: {search_term}')
    # print(f'Method: {method}')
    if tree is not None:
        tree.destroy()
    table_frame = ttk.Frame(second_tab)
    table_frame.grid(row=5, column=0, padx=(10,0), sticky=tk.NSEW) 
    table_frame.grid_rowconfigure(0, weight=1)  # Add this line
    table_frame.grid_columnconfigure(0, weight=1)  # Add this line
    tree = ttk.Treeview(table_frame, height=10)
    
    

    if method == 'Exact Match':
        # try:
            exact_matches = exact_match(search_term, substrate_df, termini_list, output_directory)
            matches_count = len(exact_matches)
            # Clear previous data in the treeview
            for i in tree.get_children():
                tree.delete(i)
            columns = list(exact_matches.columns)
        # Create the columns in the tree widget.
            tree["columns"] = columns
            for column in columns:
                tree.column(column, stretch=False, anchor=tk.CENTER, width=125)
                tree.heading(column, text=column)
            # Add the data to the tree table
            for index, row in exact_matches.iterrows():
                tree.insert('', 'end', values=list(row))
            # return exact_matches
        # except Exception as e:
        #     messagebox.showerror("Error", f"An unexpected error occurred while attempting to match termini to cleavage sites.")
        #     return
    

    elif method == 'Fuzzy Match':
        # try:
            fuzzy_matches = fuzzy_match(search_term, substrate_df, termini_list, output_directory)
            matches_count = len(fuzzy_matches)
            for i in tree.get_children():
                tree.delete(i)
            columns = list(fuzzy_matches.columns)
        # Create the columns in the tree widget.
            tree["columns"] = columns
            for column in columns:
                tree.column(column, stretch=False, anchor=tk.CENTER, width=125)
                tree.heading(column, text=column)
            # Add the DataFrame data to the tree widget.
            for index, row in fuzzy_matches.iterrows():
                tree.insert('', 'end', values=list(row))
        #     return fuzzy_matches
        # except Exception as e:
        #     messagebox.showerror("Error", f"An unexpected error occurred while attempting to match termini to cleavage sites.")
        #     return
    

    scrollbar_x = ttk.Scrollbar(table_frame, orient=tk.HORIZONTAL)
    scrollbar_x.config(command=tree.xview)
    scrollbar_x.grid(row=6, sticky=tk.EW)
    tree.config(xscrollcommand=scrollbar_x.set)

    scrollbar_y = ttk.Scrollbar(table_frame, orient=tk.VERTICAL)
    scrollbar_y.config(command=tree.yview)
    scrollbar_y.grid(row=5, column=2, sticky=tk.NS)
    tree.config(yscrollcommand=scrollbar_y.set) 
    tree.config(selectmode="extended")
    tree.grid(row=5, sticky=tk.NSEW)
    
    # second_tab.grid_rowconfigure(5, weight=1)
    # second_tab.grid_columnconfigure(0, weight=1)
    tree['show'] = 'headings'
    table_label["text"] = f"Found {matches_count} matches using {method}."

    # def on_resize(event):
    #     # Calculate the total treeview width
    #     # width = event.width
    #     width = second_tab.winfo_width()
    #     n = len(tree["columns"])   # include identifier column as well
    #     # Set each column's width to the total width divided by the number of columns
    #     for column in tree["columns"]:
    #         tree.column(column, width=width//n)
    # second_tab.bind('<Configure>', on_resize)







    


######################
# Window and widgets #
######################
window = tk.Tk()
window.title("Non-Tryptic Peptide Extractor")
window.geometry("1000x700")
window.resizable(width=True, height=False)


# Style configuration of widgets
style = ttk.Style()
window.tk.call("source", "Azure-ttk-theme-main/azure.tcl")
window.tk.call("set_theme", "dark")

style.configure("TButton", font=('TkDefaultFont'))
style.configure("TListbox", padding=10)
style.configure("TProgressBar", thickness=50)
style.layout("TNotebook", [])
style.configure("TNotebook", tabmargins=0)


notebook = ttk.Notebook(window, style="TNotebook")
notebook.grid()
main_tab = ttk.Frame(notebook, width=1200, height=1600)
main_tab.grid()

notebook.add(main_tab, text="Non-Tryptic Peptide Extraction")

# scrollbar = tk.Scrollbar(window)
# scrollbar.grid()

# Define + arrange widgets in the interface

#### Main tab ####
frame_browse = ttk.Frame(main_tab)
frame_browse.grid(row=2, column=0, sticky=tk.W)
browse_button = ttk.Button(frame_browse, text="Browse Files", command=browse_files, style='Accent.TButton')
browse_button.grid(row=2, column=0, pady=10, padx=10, sticky=tk.W)

frame_remove = ttk.Frame(main_tab)
frame_remove.grid(row=2, column=1, sticky=tk.W)
remove_button = ttk.Button(frame_browse, text="Remove Selected", command=remove_files, style='Accent.TButton')
remove_button.grid(row=2, column=1, pady=10, padx=10, sticky=tk.W)

frame_process = ttk.Frame(main_tab)
process_button = ttk.Button(main_tab, text="Process", command=process_files_button, state=tk.NORMAL)
process_button.grid(row=4, column=0, pady=10, padx=10, sticky=tk.W)

listbox_label = tk.Label(main_tab, text="List of selected PSM files:")
listbox_label.grid(row=0, column=0, pady=10, padx=10, sticky=tk.W)

listbox_frame = tk.Frame(main_tab)
listbox_frame.grid(row=1, column=0, columnspan=2, sticky=tk.W)
# file_listbox = tk.Listbox(listbox_frame, selectmode=tk.MULTIPLE, width=100, height=7, highlightbackground="#008BFF")
# file_listbox.grid(row=1, column=0, columnspan=2, padx=(10,0), sticky=tk.W)
# scrollbar = ttk.Scrollbar(listbox_frame, orient="vertical")
# scrollbar.config(command=file_listbox.yview)
# scrollbar.grid(row=1, column=2, sticky=tk.NS)
# file_listbox.config(yscrollcommand=scrollbar.set)

completion_label = tk.Label(main_tab, text="", fg='#007fff', font=("TkDefaultFont", 11, "bold"), justify= tk.LEFT)
completion_label.grid(row=9, column=0, pady=10, padx=5, sticky=tk.W)

frame_output_files = tk.Frame(main_tab)
frame_output_files.grid(row=12, column=0, sticky=tk.W)
output_files_label = tk.Label(main_tab, text="", justify= tk.LEFT)
output_files_label.grid(row=12, column=0,pady=10, padx=10, sticky='w')

frame_output_directory = tk.Frame(main_tab)
frame_output_directory.grid(row=3, column=0, sticky=tk.W)
output_directory_button = ttk.Button(frame_output_directory, text="Select Output Directory", command=select_output_directory, style='Accent.TButton')
output_directory_button.grid(row=3, column=0, pady=10, padx=10, sticky=tk.W)

frame_output_directory_label = tk.Frame(main_tab)
frame_output_directory_label.grid(row=3, column=1, sticky=tk.W)
output_directory_label = tk.Label(frame_output_directory, text="No output directory selected yet.")
output_directory_label.grid(row=3, column=1, pady=10, padx=10, sticky=tk.W)

barnum = tk.IntVar()
loading_bar = ttk.Progressbar(main_tab, orient='horizontal', mode='determinate', length=200, variable=barnum, maximum=100.0)
loading_bar.grid(row=8, column=0, pady=10, padx=10, sticky=tk.W)

group_button = ttk.Button(frame_browse, text="Group Files", command=group_files, style='Accent.TButton')
group_button.grid(row=2, column=2, pady=10, padx=10, sticky=tk.W)

# group_mutant_button = ttk.Button(main_tab, text="Group as Mutant", command=group_as_mutant, style='Accent.TButton')
# group_mutant_button.grid(row=4, column=1, pady=10, padx=10, sticky=tk.W)

# group_listbox = tk.Listbox(listbox_frame, selectmode=tk.MULTIPLE, width=100, height=7)
# group_listbox.grid(row=2, column=0, columnspan=2, padx=(10,0), pady=10, sticky=tk.W)
# Initialize the treeview
treeview = ttk.Treeview(listbox_frame, columns=("File Path", "Group"), show="headings", height=5)
treeview.heading("File Path", text="File Path")
treeview.heading("Group", text="Group")
treeview.column("File Path", width=600)
treeview.column("Group", width=75, anchor=tk.CENTER)

treeview.grid(row=1, column=0, padx=(10,0), sticky=tk.W)

# Add scrollbar
scrollbar = ttk.Scrollbar(listbox_frame, orient="vertical")
scrollbar.config(command=treeview.yview)
scrollbar.grid(row=1, column=2, sticky=tk.NS)
treeview.config(yscrollcommand=scrollbar.set)

# scrollbar = ttk.Scrollbar(listbox_frame, orient="vertical")
# scrollbar.config(command=file_listbox.yview)
# scrollbar.grid(row=1, column=2, sticky=tk.NS)


###############################
#### Protease Matching Tab ####
second_tab = tk.Frame(notebook, width=1200, height=1600, highlightthickness=0)
second_tab.grid()
notebook.add(second_tab, text="Protease Matching")

info_label_frame = ttk.Frame(second_tab)
info_label_frame.grid(row=1, column=0, sticky=tk.W)
info_label = tk.Label(info_label_frame, text="Please perform non-tryptic peptide extraction before protease matching.", justify= tk.LEFT)
info_label.grid(row=1, column=0, pady=10, padx=10, sticky=tk.W)

menu_frame = tk.Frame(second_tab)
menu_frame.grid(row=2, column=0, sticky=tk.W)

organism_label = tk.Label(menu_frame, text="Select species:", justify= tk.LEFT)
organism_label.grid(row=2, column=0, padx=10, sticky=tk.W)
# threading.Thread(target=substrate_table_processing).start()

method_label = tk.Label(menu_frame, text="Select matching method:", justify= tk.LEFT)
method_label.grid(row=2, column=1, padx=10, sticky=tk.W)

organisms_var = tk.StringVar()
organisms_var.trace("w", update_organism_menu)
organism_menu = ttk.Combobox(menu_frame, textvariable=organisms_var, width=30, justify= tk.LEFT, values = organisms, state=tk.DISABLED)
organism_menu.grid(row=3, column=0, padx=10, pady=10, sticky=tk.W)
organism_menu.bind('<<ComboboxSelected>>', update_organism_menu)

protease_button = ttk.Button(second_tab, text="Match Termini", command=protease_match_button, state=tk.DISABLED)
protease_button.grid(row=4, column=0, padx=10, pady=10, sticky=tk.W)

method_var = tk.StringVar()
method_var.trace("w", update_method_menu)
match_method_menu = ttk.Combobox(menu_frame, textvariable=method_var, width=30, justify= tk.LEFT, values=['Exact Match', 'Fuzzy Match'], state=tk.DISABLED)
match_method_menu.grid(row=3, column=1, padx=10, pady=10, sticky=tk.W)
match_method_menu.bind('<<ComboboxSelected>>', update_method_menu)

table_label = tk.Label(second_tab, text = "", fg='#007fff', font=("TkDefaultFont", 11, "bold"), justify=tk.LEFT)
table_label.grid(row=6, column=0, padx=10, sticky=tk.W)


window.mainloop()


# if __name__ == '__main__':
#     main()

###############################

# frame_output_files = tk.Frame(window)
# frame_output_files.grid(row=6, column=0, sticky=tk.S)
# output_listbox_label = tk.Label(frame_output_files, text="Output files:")
# output_listbox_label.grid(row=6, column=0, pady=10, padx=10, sticky=tk.S)

# frame_output_files = tk.Frame(window)
# frame_output_files.grid(row=10, column=0, sticky=tk.W)
# output_files_label = tk.Label(frame_output_files, text="", anchor="w")
# output_files_label.grid(row=10, column=0, pady=10, padx=10, sticky=tk.W+tk.E)

# output_listbox = tk.Listbox(.window, selectmode=tk.MULTIPLE, width=100, height=10)
# output_listbox.grid(row=6, column=0, columnspan=2, padx=10, sticky=tk.W)