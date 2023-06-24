import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt

peptide_df_all = pd.read_csv("./peptide_output.csv")

# subset of peptide_df that only contains peptides that appear in at least 3 of the 4 mutant samples, but don't appear in WT samples
subset_df = peptide_df_all[(peptide_df_all.iloc[:, 7:11].sum(axis=1)  >= 3) & (peptide_df_all.iloc[:, 11:].sum(axis=1) == 0)]
subset_df.to_csv('./mutant_only_subset.csv', index=False)

# subset the non-tryptic peptides from peptide_df_all and from subset_df
peptide_df_nontryp = peptide_df_all[peptide_df_all['Tryptic State'] == 'Non-Tryptic']
print(len(peptide_df_nontryp)) # 312 fully non-tryptic peptides

subset_df_nontryp = subset_df[subset_df['Tryptic State'] == 'Non-Tryptic']
print(len(subset_df_nontryp)) # 6 fully non-tryptic peptides in subset df


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
# Hydrophobic AAs are A, V, I, L, P, F, G, M

# What do we want:
# Counts for Prev AA, Next AA, First AA, Last AA

def get_first_letter(text):
    return text[0]

def get_last_letter(text):
    return text[-1]


def aa_counts(peptide_df, name):
    amino_acids_dict = {
        'A': '#F5B107',
        'R': '#0D50C5',
        'N': '#0D50C5',
        'D': '#0D50C5',
        'C': '#0D50C5',
        'E': '#0D50C5',
        'Q': '#0D50C5',
        'G': '#0D50C5',
        'H': '#0D50C5',
        'I': '#F5B107',
        'L': '#F5B107',
        'K': '#0D50C5',
        'M': '#F5B107',
        'F': '#F5B107',
        'P': '#F5B107',
        'S': '#0D50C5',
        'T': '#0D50C5',
        'W': '#F5B107',
        'Y': '#0D50C5',
        'V': '#F5B107',
        '-': 'red'
    }

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

    fig, axes = plt.subplots(2,2, figsize=(9,9))

    for ax_idx, col in zip(axes.ravel(), counts_df.columns):
        counts_df_sorted = counts_df.sort_values(by=col, ascending=False)
        colors = [amino_acids_dict[amino] for amino in counts_df_sorted.index]

        # ax_idx = axes[ax]
        ax_idx.bar(counts_df_sorted.index, counts_df_sorted[col], color=colors)
        ax_idx.set_title(col)
        ax_idx.set_xlabel('Amino Acid')
        ax_idx.set_ylabel('Frequency')

        # for i, v in enumerate(col):
        #     ax_idx.text(i, v, str(v), ha='center', va='bottom')

    plt.subplots_adjust(hspace=0.5, wspace=0.5)
    plt.savefig(f'./{name}_counts_bargraph.png') 
    plt.show()

aa_counts(peptide_df_all, "all_peptides")
aa_counts(subset_df, "mutant_only_subset")