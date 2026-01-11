#%%

import os
import requests

# Input / output
proteome_id_file = "../Data/proteome_ids.txt"
output_fasta = "../Data/all_fastas.fasta"

base_url = "https://rest.uniprot.org/uniprotkb/stream"

# Retrieve and combine FASTA files
with open(output_fasta, "w") as combined_fasta:
    with open(proteome_id_file) as f:
        for line in f:
            proteome_id = line.strip()
            if not proteome_id:
                continue

            params = {
                "format": "fasta",
                "query": f"(proteome:{proteome_id})"
            }

            response = requests.get(base_url, params=params, timeout=60)

            if response.status_code == 200 and response.text.strip():
                combined_fasta.write(response.text.rstrip() + "\n")
                print(f"Added {proteome_id}")
            else:
                print(f"Failed to download {proteome_id} (status {response.status_code})")

# %%
