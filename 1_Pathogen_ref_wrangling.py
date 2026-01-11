#%%
import pandas as pd
import numpy as np
    
################ Initial filtering and wrangling ###############################

# Pathogen reference data
pathogen_ref = pd.read_csv("../Data/pathogen_ref.tsv", sep='\t')

# Proteome data
proteome = pd.read_csv("../Data/referenced_proteomes.tsv", sep='\t')

# Remove NA's
proteome = proteome.dropna(subset=["Proteome Id"])
proteome = proteome.dropna(subset=["Genome assembly ID"])

################ Specific filtering and wrangling ###############################

# Pattern for finding human pathogens
pattern = r'(?i)human|homo'

# Remove NAs
pathogen_ref = pathogen_ref.dropna(subset=['Host'])
pathogen_ref = pathogen_ref.dropna(subset=['Assembly'])

# Remove pathogens that are not associated with humans
filtered_df = pathogen_ref[pathogen_ref['Host'].str.contains(pattern, regex=True)]

# Remove the "non-human" pathogens since these escape pattern
filtered_pathogen_ref = filtered_df[~filtered_df['Host'].str.contains(r'(?i)non', regex=True)]

# Merging proteome data onto pathogen ref
merged_proteome = proteome.merge(filtered_pathogen_ref, how="left", left_on="Genome assembly ID", right_on="Assembly", suffixes=('', '_ref'))

# Remove weird instances with missing Proteome id
# Looks at NAs
unique_merged = merged_proteome.dropna(subset=['Assembly'])

# Saving certain columns from merged df as .csv
unique_merged[["Proteome Id", "#Organism group", "Strain", "Taxonomic lineage", "Protein count", "Assembly"]].to_csv("../Data/Pathogenic_bacteria_proteome.csv", index=False)


# Saving proteomes IDs as .txt for download of full proteomes later
unique_merged['Proteome Id'].to_csv("../Data/proteome_ids.txt", header=False, index=False)

# Finding annotated total protein count
print(unique_merged.shape)
print(unique_merged["Protein count"].sum())

# Check how many IDs match before merging (number of human pathogenic proteomes)
matching_ids = set(proteome["Genome assembly ID"]) & set(filtered_pathogen_ref["Assembly"])
print(f"Number of matching Assembly IDs: {len(matching_ids)}")


# %%
