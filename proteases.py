import pandas as pd
import os
import numpy as np
import time
import re
from bioservices import UniProt
import requests
from unipressed import IdMappingClient

all_proteins_df = pd.read_csv('./all_unique_proteins.csv')
proteases_df = all_proteins_df[all_proteins_df['Protein Description'].str.contains(r'(?:protease|peptidase)')]
    
proteases_df.to_csv('./proteases_df.csv', index=False)
# Read in file containing all unique proteins across samples
def protease_func(df):

    # subset proteases/peptidases
    # proteases_df = all_proteins_df[all_proteins_df['Protein Description'].str.contains(r'(?:protease|peptidase)')]

    # proteases_df.to_csv('./proteases_df.csv', index=False)

    # create csv file that contains accession IDs for the proteases that will then be used to map to the MEROPS database.
    output_file = 'accession_IDs.csv'
    merops_id_list = []
    # Create list of accessions
    accessions = [str(value) for value in df['Protein ID']]

    # Iterate through list of accession IDs
    for accession in accessions:
        url = f"https://www.uniprot.org/uniprot/{accession}.txt"
        # Send HTTP GET request to the UniProt url to retrieve content of entry for accession ID
        response = requests.get(url)
        # Checks if response from server is successful before going forward
        if response.ok:
            lines = response.text.split("\n")
            # look for line that starts with DR MEROPS that contains the MEROPS ID for this protease. If matching line found, that line is assigned to the mapping_line variable, otherwise it equals None
            mapping_line = next((line for line in lines if line.startswith("DR   MEROPS;")), None)
            if mapping_line:
                # regex to match MEROPS; in the line
                mapping = re.search(r"MEROPS;\s*(\S+?)(?:;|$)", mapping_line)
                if mapping:
                    merops_id1 = mapping.group(1)
                    merops_id_list.append(merops_id1)
                    print(f"Protein ID: {accession}, MEROPS ID: {merops_id1}")
                else:
                    merops_id_list.append(None)
            else:
                    merops_id_list.append(None)
            time.sleep(1)
    print(merops_id_list)
    df.loc[:,'MEROPS ID'] = merops_id_list
    print(df)
    # merops_id = ''.join(merops_id_list)
    # return pd.Series({'MEROPS ID': merops_id})
    return df



   
    # for id in accessions:
    #     request = IdMappingClient.submit(
    #         source="UniProtKB_AC-ID", dest="MEROPS", ids=id
    #     )
    #     time.sleep(5)
    #     print(list(request.each_result()))
protease_func(proteases_df)
# proteases_df[['MEROPS ID']] = proteases_df.apply(protease_func, axis=1)
# proteases_df.to_csv('./protease_df_merops.csv', index=False)

# acc = ['W1HUC5','W1HT22','W1HLM8','W1HPD9','W1HX15', 'W1I2E8']

 # with open(output_file, 'w') as file:
    #     file.write(accessions)
    # print(accessions)