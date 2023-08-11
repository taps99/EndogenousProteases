# Non-Tryptic Peptide Extractor and Protease Mapping

GUI developed using Tkinter framework in Python to extract missed cleavage data from processed proteomics data resulting from Fragpipe and match termini to cleavage sites in the MEROPS database.




### Output files from non-tryptic peptide extraction
Row counts are based on the validation data set used for the development of this project (*Klebsiella pneumoniae* wild-type vs Lon mutant data from Muselius, B., Sukumaran, A., Yeung, J., & Geddes-McAlister, J. (2020). Iron Limitation in Klebsiella pneumoniae Defines New Roles for Lon Protease in Homeostasis and Degradation by Quantitative Proteomics. Frontiers in Microbiology, 11. https://www.frontiersin.org/articles/10.3389/fmicb.2020.00546)

* 1. all_peptides.csv (~301,379 rows)
- Contains all peptides across all samples. 

2. all_unique_proteins.csv (1974 rows)
- Contains all unique proteins across all samples (tryptic + non-tryptic) based on accession ID

3. all_proteases.csv
- Contains all unique proteases across all samples (unique based on Protein ID, proteins that contain "protease" or "peptidase" in their names) along with presence/absence (P/A) data

4. non-tryptic_unique_proteins.csv (672 rows)
- Contains P/A data the unique proteins associated with non-tryptic peptides found across all samples.

5. unique_non-tryptic_sequences (1331 rows)
- Contains all unique non-tryptic peptides along with the protein along with protein sequence from UniProt for associated protein. Also contains columns containing non-tryptic termini.

6. non_unique_non-tryptic_sequences (3438 rows)
- Contains ALL non-tryptic peptides (w/ duplicates) along with protein sequence from UniProt for associated protein. Also contains columns containing non-tryptic termini.

7. unique_non-tryptic_termini
- Contains a list of unique non-tryptic termini generated from non-tryptic peptides and their associated protein sequences. Also contains P/A data.

8. non_unique_non-tryptic_termini
- Contains a list of all non-tryptic termini generated from non-tryptic peptides and their associated protein sequences. Also contains P/A data.