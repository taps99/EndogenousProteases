import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
import inspect

peptide_df_all = pd.read_csv("./peptide_output.csv")

# subset of peptide_df that only contains peptides that appear in at least 3 of the 4 mutant samples, but don't appear in WT samples
subset_df_4 = peptide_df_all[(peptide_df_all.iloc[:, 7:11].sum(axis=1)  > 3) & (peptide_df_all.iloc[:, 11:].sum(axis=1) == 0)]

# Output accession IDs to file
    # output_file = 'accession_IDs.csv'
    # accessions = ','.join(str(value) for value in peptide_df['Protein ID'])
    # with open(output_file, 'w') as file:
    #     file.write(accessions)

# print(len(peptide_df['Protein'].unique())) # 632 unique proteins, 672 unique protein IDs?
# print(len(subset_df_4['Protein'].unique()))
# print((subset_df_4['Protein'].unique())) # 23 unique proteins that are only present in mutant
# counts = subset_df_4['Prev AA'].value_counts()


# The hydrophobic amino acids include alanine (Ala, A), valine (Val, V), leucine (Leu, L), isoleucine (Ile, I), proline (Pro, P), phenylalanine (Phe, F) and cysteine (Cys, C).
# Hydrophobic AAs are A, V, L, I, P, F, C

# What do we want:
# Counts for Prev AA, Next AA, First AA, Last AA

def get_first_letter(text):
    return text[0]

def get_last_letter(text):
    return text[-1]


def aa_counts(peptide_df, name):
    peptide_df['First AA'] = peptide_df['Peptide'].apply(get_first_letter)
    peptide_df['Last AA'] = peptide_df['Peptide'].apply(get_last_letter)


    prev_AA_counts = peptide_df['Prev AA'].value_counts()
    next_AA_counts = peptide_df['Next AA'].value_counts()
    first_AA_counts = peptide_df['First AA'].value_counts()
    last_AA_counts = peptide_df['Last AA'].value_counts()

    # peptide_df.to_csv("peptide_subset_Mutant_only.csv", index=False)

    counts_df = pd.DataFrame({'Prev AA': prev_AA_counts, 'First AA': first_AA_counts, 'Next AA':next_AA_counts, 'Last AA':last_AA_counts})
    counts_df.index.name = 'Amino Acid'
    counts_df.to_csv(f'./{name}_counts_df.csv')

    amino_acids_dict = {
        'A': 'hydrophobic',
        'R': 'hydrophilic',
        'N': 'hydrophilic',
        'D': 'hydrophilic',
        'C': 'hydrophobic',
        'E': 'hydrophilic',
        'Q': 'hydrophilic',
        'G': 'hydrophobic',
        'H': 'hydrophilic',
        'I': 'hydrophobic',
        'L': 'hydrophobic',
        'K': 'hydrophilic',
        'M': 'hydrophobic',
        'F': 'hydrophobic',
        'P': 'hydrophobic',
        'S': 'hydrophilic',
        'T': 'hydrophilic',
        'W': 'hydrophobic',
        'Y': 'hydrophilic',
        'V': 'hydrophobic'
    }

    color_dict = {
        'hydrophobic': '#F5B107',
        'hydrophilic': '#0D50C5'
    }

    colors = [color_dict.get(amino_acids_dict.get(aa, 'hydrophilic')) for aa in counts_df.index]

    fig, axes = plt.subplots(2,2, figsize=(9,9))

    for ax_idx, col in zip(axes.ravel(), counts_df.columns):
        counts_df_sorted = counts_df.sort_values(by=col, ascending=False)
        # ax_idx = axes[ax]
        ax_idx.bar(counts_df_sorted.index, counts_df_sorted[col], color=colors)
        ax_idx.set_title(col)
        ax_idx.set_xlabel('Amino Acid')
        ax_idx.set_ylabel('Frequency')

        # for i, v in enumerate(col):
        #     ax_idx.text(i, v, str(v), ha='center', va='bottom')

    # # Create the bar chart using Pandas plotting
    # prev_AA_counts.plot(kind='bar', ax=ax1)
    # ax1.set_title('Prev AA')
    # ax1.set_xlabel('Amino Acid')
    # ax1.set_ylabel('Frequency')

    # for i, v in enumerate(prev_AA_counts):
    #     ax1.text(i, v, str(v), ha='center', va='bottom')

    # first_AA_counts.plot(kind='bar', ax=ax2)
    # ax2.set_title('First AA')
    # ax2.set_xlabel('Amino Acid')
    # ax2.set_ylabel('Frequency')

    # for i, v in enumerate(first_AA_counts):
    #     ax2.text(i, v, str(v), ha='center', va='bottom')

    plt.subplots_adjust(hspace=0.5, wspace=0.5)
    plt.savefig(f'./{name}_counts_bargraph.png') 


    # Display the chart
    plt.show()

aa_counts(peptide_df_all, "all_peptides")
aa_counts(subset_df_4, "mutant_only_subset")