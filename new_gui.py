import tkinter as tk
from tkinter import ttk
from tkinter import filedialog
from tkinter import messagebox
from tkinter import simpledialog
import pandas as pd
import numpy as np
import os
import sys
from new_functions import process_data, peptide_seq_match, get_protein_sequence, output_files, create_termini_list
import threading
import time
from matching import exact_match, fuzzy_match

def main():
    global substrate_df
    global output_directory
    # Check if the application is running as a script or a packaged application
    if getattr(sys, 'frozen', False):
    # The application is running as a bundled executable
        application_path = sys._MEIPASS
    else:
    # The application is running as a script
        application_path = os.path.dirname(os.path.abspath(__file__))
    # Set up path variables for substrates and organism lists
    substrate_path = os.path.join(application_path, 'substrates_merops.csv')
    organisms_path = os.path.join(application_path, 'organisms_merops.csv')
    substrate_df = pd.read_csv(substrate_path)
    organisms = pd.read_csv(organisms_path)
    organisms = organisms.iloc[:,0].tolist() # convert to list for dropdown menu
    global termini_file_path
    global filtered_termini_file_path
    termini_file_path = None
    output_directory = None
    filtered_termini_file_path = None

    # Function to open file explorer to select files to be processed
    def browse_files():
        file_paths = filedialog.askopenfilenames(filetypes=[("TSV files", "*.tsv"), ("All files", "*.*")], multiple=True)
        for file_path in file_paths:
            # file_listbox.insert(tk.END, file_path)
            treeview.insert("", tk.END, values=(file_path, ""))
        switch_button_state()

    # Function to remove selected files from list
    def remove_files():
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
            output_directory_label_2["text"] = f"Output directory: {output_directory}"
        else: 
            output_directory_label["text"] = f"No output directory selected."
            output_directory_label_2["text"] = f"No output directory selected."
        os.makedirs(output_directory, exist_ok=True)
        switch_button_state()
        return output_directory

    # Initialize group names for presence/absence data
    def setup():
        count = 0
        processed_filepaths = []
        listbox_entries = []

        group_counts = {} 
        for value in treeview.get_children():
            file_path, group = treeview.item(value, "values")
            listbox_entries.append(file_path)
            file_name = os.path.basename(file_path)
            # Check if the group is already in group_counts
            if not group:  
                group = file_name
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

    # Group files button function
    def group_files():
        global group_name
        if len(treeview.selection()) > 0:
            group_name = simpledialog.askstring("Input", "Enter the group name:")
            if group_name:
                selected_items = treeview.selection()
                for item in selected_items:
                    file_path, _ = treeview.item(item, "values")
                    treeview.item(item, values=(file_path, group_name))
        else:
            messagebox.showerror("Error", "No files were selected.")


    # Process button function
    def process_files_button():
        global process_thread
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
            proceed = messagebox.askyesno("Warning", "The input files will be grouped by file name by default. Would you like to proceed?")
            if not proceed:
                return

        process_thread = threading.Thread(target=processing_functions)
        process_thread.start()

    # Switch state of process button once at least one file and output directory has been selected.
    def switch_button_state():
        if len(treeview.get_children()) > 0 and output_directory != "":
            process_button['style'] = 'Accent.TButton'
        else:
            process_button['style'] = 'TButton'

    # This function is run on another thread, which holds all the main data processing functions for non-tryptic peptide extraction
    def processing_functions():
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
            all_non_tryp_peptides = output_files(processed_filepaths, processed_dfs, unprocessed_dfs, output_directory)
        except Exception as e:
            messagebox.showerror("Error", f"An unexpected error occurred while generating the output files: {str(e)}")
            completion_label["text"] = "Process aborted. Please try again."
            barnum.set(0)
            return
        for i in range(35, 50):
            barnum.set(i)
            time.sleep(0.05)
        completion_label["text"] = f"Retrieving protein sequences from UniProt for {str(len(all_non_tryp_peptides))} peptides..."
        for i in range(50, 81):
            barnum.set(i)
            time.sleep(0.05)
        try:
            all_peptide_df, peptide_df = peptide_seq_match(get_protein_sequence(all_non_tryp_peptides), "non-tryptic_sequences", output_directory)
        except Exception as e:
            messagebox.showerror("Error", f"An unexpected error occurred while acquiring sequences from UniProt: {str(e)}")
            completion_label["text"] = "Process aborted. Please try again."
            barnum.set(0)
            return
        termini_list = create_termini_list(peptide_df,"unique_non-tryptic_termini", output_directory)
        termini_list_all = create_termini_list(all_peptide_df, "non_unique_non-tryptic_termini", output_directory)
        completion_label["text"] = "Finalizing..."
        for i in range(80, 101):
            barnum.set(i)
            time.sleep(0.05)

        count = len(listbox_entries)
        if count == 1:
            completion_label["text"] = f"Successfully processed {count} file."
        else:
            completion_label["text"] = f"Successfully processed {count} files."

        output_filepaths = [f'{output_directory}/all_peptides.csv', f'{output_directory}/all_unique_proteins.csv', f'{output_directory}/all_proteases.csv', 
                            f'{output_directory}/non-tryptic_unique_proteins.csv', f'{output_directory}/unique_non-tryptic_sequences.csv', f'{output_directory}/non_unique_non-tryptic_sequences.csv', f'{output_directory}/unique_non-tryptic_termini.csv', f'{output_directory}/non_unique_non-tryptic_termini.csv']
        # output_filepaths = [os.path.join(output_directory,f) for (dirpath, dirnames, filenames) in os.walk(output_directory) for f in filenames]
        # for dirpath, dirnames, filenames in os.walk(output_directory):
        #     for filename in filenames:
        #         output_filepaths.append(os.path.join(dirpath, filename))
        paths_text = '\n'.join(output_filepaths)
        output_files_label["text"] = f"Output files:\n{paths_text}"

        return termini_list, termini_list_all


    def select_termini_file():
        global termini_file_path
        global termini_count
        termini_file_path = filedialog.askopenfilename(filetypes=[("CSV files", "*.csv"),("All files", "*.*")], multiple=False)
        df = pd.read_csv(termini_file_path)
        termini_count = str(len(df))
        if len(str(termini_file_path)) > 0:
            termini_file_label['text'] = f'File: {termini_file_path}'
            termini_count_label["text"] = f'Matching with {termini_count} termini.'
        else:
            termini_file_label['text'] = "No file has been selected. Please select a file containing non-tryptic termini."
        return termini_file_path

    def open_groupfilter_window():
        global termini_file_path
        global output_directory
        # global selected_option_dict
        # global group_name_list
        # if termini_file_path:
        if not output_directory:
            messagebox.showerror("Error", "Please select output directory before grouping.")
            return
        try:
            df = pd.read_csv(termini_file_path) 
        except Exception as e:
            messagebox.showerror("Error", "Please input a .csv file containing non-tryptic termini.")
            return
        if len(list(df.columns)) <= 1:
            messagebox.showerror("Error", "File does not contain presence/absence data.")
            return

        groups = [] # list of tuples (group_name, column) -> (MUTANT, MUTANT_1), (MUTANT, MUTANT_2), etc.
        for col in df.columns[4:]:
            group_name = '_'.join(col.rsplit('_')[:-1])
            groups.append((group_name, col))

        new_window = tk.Toplevel(window)
        new_window.geometry('500x500')
        new_window.grid_columnconfigure(0, weight=1)  # Adds empty space on left
        new_window.grid_columnconfigure(3, weight=1)  # Adds empty space on right
        row_count = 1
        check_frame = tk.Frame(new_window)
        check_frame.grid(sticky=tk.NSEW)

        # group_label = ttk.Label(check_frame, text="Groups:")
        # group_label.grid(row=0)

        def create_checkbutton(var, dropdown1, dropdown2):
            return lambda: toggle_dropdown(var, dropdown1, dropdown2)
        # grab first element of each tuple in the groups list and put them into a list
        unique_groups = list(set([x[0] for x in groups]))
        selected_option_dict = {}
        
        for group in unique_groups: # iterate through each unique group 
            
            blah = tk.IntVar()
            options = [0] + [i for i in range(0, sum([1 for x in groups if x[0] == group]) + 1)]
            selected_option = tk.IntVar()
            selected_option_dict[group] = selected_option

            dropdown = ttk.OptionMenu(new_window, selected_option, *options)
            dropdown.configure(state='disabled')
            dropdown.grid(row=row_count, column=3, pady=(10,0), padx=(0,10), sticky=tk.W)

            comparator_options = ['Equal to', 'Less Than or Equal to', 'More Than or Equal to']
            selected_comparator = tk.StringVar()  # create a StringVar for the comparator
            selected_option_dict[group] = (selected_option, selected_comparator)  # store tuple of selected option and comparator

            comparator_dropdown = ttk.OptionMenu(new_window, selected_comparator, comparator_options[0], *comparator_options)
            comparator_dropdown.configure(state='disabled')
            comparator_dropdown.grid(row=row_count, column=2, pady=(10,0), sticky=tk.W)  # adjust column index as needed

            check = ttk.Checkbutton(new_window, text=group, variable=blah, onvalue=1, offvalue=0, command=create_checkbutton(blah, dropdown, comparator_dropdown))
            check.grid(row=row_count, column=1, padx=(0,20), pady=(10,0), sticky=tk.W)
            row_count += 1  
        
        def toggle_dropdown(var, dropdown1, dropdown2):
            if var.get() == 1: # if the checkbutton is selected
                dropdown1.configure(state='normal')
                dropdown2.configure(state='normal')
            else: # if the checkbutton is deselected
                dropdown1.configure(state='disabled')
                dropdown2.configure(state='normal')


        # Function to apply the filter when the "Apply filter" button is clicked
        def apply_filter(): 
            # global selected_option_dict
            global termini_file_path
            global filtered_termini_file_path
            df = pd.read_csv(termini_file_path)
            masks = []
            for group, (option_var,comparator_var) in selected_option_dict.items():  
                # find all columns belonging to each group 
                columns_for_group = [col for group_name, col in groups if group_name == group]
                # then retain only rows where sum of the values in group's columns is >= to the number selected in dropdown menu for that group
                # mask = df[columns_for_group].sum(axis=1) >= option_var.get()
                if comparator_var.get() == 'Equal to':
                    mask = df[columns_for_group].sum(axis=1) == option_var.get()
                elif comparator_var.get() == 'Less Than or Equal to':
                    mask = df[columns_for_group].sum(axis=1) <= option_var.get()
                elif comparator_var.get() == 'More Than or Equal to':
                    mask = df[columns_for_group].sum(axis=1) >= option_var.get()
                masks.append(mask)
            
            final_mask = np.logical_and.reduce(masks)
            df = df[final_mask]
            print(len(df))
            termini_count = str(len(df))
            termini_count_label["text"] = f"Matching with {termini_count} termini."
            df.to_csv(f'{output_directory}/filtered_termini.csv', index=False)
            filtered_termini_file_path = f'{output_directory}/filtered_termini.csv'
            # termini_file_path = filtered_termini_file_path
            new_window.destroy()
            return termini_file_path
        apply_frame = ttk.Frame(new_window)
        apply_frame.grid(row=row_count+1)
        apply_button = ttk.Button(new_window, text="Apply Group Filter", command=lambda: apply_filter(), style="Accent.TButton")
        apply_button.grid(row=row_count+1, column=2, pady=(30,0), padx=(0,10), sticky=tk.W)
            
        
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


    # tree = None
    def protease_match():
        global substrate_df
        global termini_list
        global tree
        tree = None
        search_term = update_organism_menu().capitalize()
        method = update_method_menu()
        if tree is not None:
            tree.destroy()
        table_frame = ttk.Frame(second_tab)
        table_frame.grid(row=5, column=0, padx=(10,0), sticky=tk.NSEW) 
        table_frame.grid_rowconfigure(0, weight=1)  # Add this line
        table_frame.grid_columnconfigure(0, weight=1)  # Add this line
        tree = ttk.Treeview(table_frame, height=10)
        

        if method == 'Exact Match':
            try:
                if filtered_termini_file_path is not None:
                    exact_matches = exact_match(search_term, substrate_df, filtered_termini_file_path, output_directory)
                else:
                    exact_matches = exact_match(search_term, substrate_df, termini_file_path, output_directory)
            except KeyError as k:
                messagebox.showerror("Error", f"This file does not contain the correct column names: {str(k)}")
                return
            except Exception as e:
                messagebox.showerror("Error", f"An unexpected error occurred while attempting to match termini to cleavage sites: {str(e)}")
                return
            matches_count = len(exact_matches)
            # Clear previous data in the treeview
            for i in tree.get_children():
                tree.delete(i)
            columns = list(exact_matches.columns)
        # Create the columns in the tree widget
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
            try:
                if filtered_termini_file_path is not None:
                    fuzzy_matches = fuzzy_match(search_term, substrate_df, filtered_termini_file_path, output_directory)
                else:
                    fuzzy_matches = fuzzy_match(search_term, substrate_df, termini_file_path, output_directory)
            except KeyError as k:
                messagebox.showerror("Error", f"This file does not contain the correct column names: {str(k)}")
                return
            except Exception as e:
                messagebox.showerror("Error", f"An unexpected error occurred while attempting to match termini to cleavage sites: {str(e)}")
                return
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

    window.after(0, update_organism_menu, organisms)
    # Style configuration of widgets
    style = ttk.Style()
    window.tk.call("source", "Azure-ttk-theme-main/azure.tcl")
    window.tk.call("set_theme", "dark")

    style.configure("TButton", font=('TkDefaultFont'))
    style.configure("TListbox", padding=10)
    style.configure("TProgressBar", thickness=50)
    style.layout("TNotebook", [])
    style.configure("TNotebook", tabmargins=0)
    style.configure("TCheckbutton")


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

    # grouping_button1 = ttk.Button(main_tab, text="Grouping Filter", command=open_groupfilter_window, state=tk.NORMAL)
    # grouping_button1.grid(row=4, column=1, sticky=tk.W)

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
    info_label_frame.grid(row=0, column=0, sticky=tk.W)
    info_label = tk.Label(info_label_frame, text="Select a file containing non-tryptic termini.", justify= tk.LEFT)
    info_label.grid(row=0, column=0, pady=10, padx=10, sticky=tk.W)

    termini_file_button = ttk.Button(info_label_frame, text="Select File", command=select_termini_file, state=tk.NORMAL, style = 'Accent.TButton')
    termini_file_button.grid(row=1, column=0, padx=10, pady=(0,10), sticky=tk.W)

    output_directory_button_2 = ttk.Button(info_label_frame, text="Select Output Directory", command=select_output_directory, style='Accent.TButton')
    output_directory_button_2.grid(row=2, column=0, pady=10, padx=10, sticky=tk.W)

    termini_file_label = ttk.Label(info_label_frame, text="", justify= tk.LEFT)
    termini_file_label.grid(row=1, column=0, padx=(110,0), pady=(0,10), sticky=tk.W)

    output_directory_label_2 = tk.Label(info_label_frame, text="No output directory selected.")
    output_directory_label_2.grid(row=2, column=0, padx=(170,0), sticky=tk.W)

    menu_frame = tk.Frame(second_tab)
    menu_frame.grid(row=2, column=0, sticky=tk.W)

    grouping_filter_button = ttk.Button(menu_frame, text="Grouping Filter", command=open_groupfilter_window, state=tk.NORMAL, style = 'Accent.TButton')
    grouping_filter_button.grid(row=3, column=2, padx=(30,0), sticky=tk.W)

    organism_label = tk.Label(menu_frame, text="Select species:", justify= tk.LEFT)
    organism_label.grid(row=2, column=0, padx=10, sticky=tk.W)
    # threading.Thread(target=substrate_table_processing).start()

    method_label = tk.Label(menu_frame, text="Select matching method:", justify= tk.LEFT)
    method_label.grid(row=2, column=1, padx=10, sticky=tk.W)

    organisms_var = tk.StringVar()
    organisms_var.trace("w", update_organism_menu)
    organism_menu = ttk.Combobox(menu_frame, textvariable=organisms_var, width=30, justify= tk.LEFT, values = organisms, state=tk.NORMAL)
    organism_menu.grid(row=3, column=0, padx=10, pady=10, sticky=tk.W)
    organism_menu.bind('<<ComboboxSelected>>', update_organism_menu)

    protease_button = ttk.Button(second_tab, text="Match Termini", command=protease_match_button, state=tk.NORMAL, style = 'Accent.TButton')
    protease_button.grid(row=4, column=0, padx=10, pady=10, sticky=tk.W)

    termini_count_label = ttk.Label(second_tab, text="", justify= tk.LEFT)
    termini_count_label.grid(row=4, column=0, pady=10, padx=(120,0), sticky=tk.W)

    method_var = tk.StringVar()
    method_var.trace("w", update_method_menu)
    match_method_menu = ttk.Combobox(menu_frame, textvariable=method_var, width=30, justify= tk.LEFT, values=['Exact Match', 'Fuzzy Match'], state=tk.NORMAL)
    match_method_menu.grid(row=3, column=1, padx=10, pady=10, sticky=tk.W)
    match_method_menu.bind('<<ComboboxSelected>>', update_method_menu)

    table_label = tk.Label(second_tab, text = "", fg='#007fff', font=("TkDefaultFont", 11, "bold"), justify=tk.LEFT)
    table_label.grid(row=6, column=0, padx=10, sticky=tk.W)


    window.mainloop()

