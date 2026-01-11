**Computationally predicting T-cell cross-reactivity induced autoimmunity with pathogenic proteome**



The purpose of the following project is to locate possible sequence

alignment between pathogenic proteins and known human autoimmune epitopes. For more in depth see attached PDF paper ("Bachelor thesis Lucio Stevns s225031.pdf").





**DATA COLLECTION:**



* Collect data from IEDB.org with epitopes and metadata of interest, this should at the very least contain protein sequences. Current script is specific for the following metadata: (Assay ID - IEDB IRI, Epitope - Name, Epitope - Molecule Parent, Epitope - Molecule Parent IRI, 1st in vivo Process - Disease, 1st in vivo Process - Disease Stage, MHC Restriction - Name, Epitope - Starting Position, Epitope - Ending Position Epitope - Modified residues). If not change script to fit metadata



* Collect data from uniprot.org with proteomes of pathogens of interest and download .tsv file. Including Assembly id is important to verify pathogenicity.

 

* Collect data from ncbi.nlm.nih.gov/pathogens/isolates (This is for verifying chosen pathogens are pathogenic). Including Assembly id is important.





**CODE STUFF**



**1\_IEDB\_Wrangling.py**



* Description: This will filter epitope of being in specific length range (12-25) and having modified residues. to reduce redundancy in epitope list the script will check if smaller epitopes exist within other longer epitopes (nested epitopes), it will do this by creating all possible 9mers for each epitope if 2 epitopes are identical it will remove the parent epitope of the 9mer that is smallest. It will then write a .csv file, with chosen columns as columns and each row being a 9mer
* Input: Autoimmune epitopes sourced from IEDB.org as .tsv file
* Output: table with unique autoimmune  



**1\_Pathogen\_ref\_wrangling**



* Description: It will clean and filter data. It will take the uniport bacterial bulk data and check how many of them have been confirmed to be pathogenic against humans, it will do this with the ncbi isolates data by cross-referencing their assembly ID.
* Input: ncbi.nlm.nih.gov/pathogens/isolates for confirmed pathogenicity and uniport.org data containing proteomes of choice.
* output: List of bacterial proteome IDs that is human pathogenic.



**2\_Pathogen\_fasta\_retrivial**



* Description: A rest api to retrive the full fasta proteomes of the bacteria.
* Input: Output from **1\_Pathogen\_ref\_wrangling (List of proteome IDs)**
* Output: 1 fasta file with all bacterial proteomes



**3\_Protein\_meta\_data**



* Description: This will in essence just get all the meta data from the bacteria included. It will be each protein in the fastafiles metadata, so which Genus-species-strain each protein comes from, protein annotation, gene name, sequence and protein\_ID.
* Input: Output from 2\_Pathogen\_fasta\_retrivial (bacterial proteomes fasta)
* Output: Table with protein as rows and metadata as columns (.csv)





**4\_Perfect\_match\_2\_0** (Weird script name i know) 



* Description: This is script does the main purpose of the project. The matching utilizes to main methods, to say it shortly it looks to see if the epitopes matches 1 to 1 with the pathogenic proteins. But since my epitopes can be up to 25 amino acids long, i utilize a sliding window approach that creates all possible sub epitopes down to 9mers, this results in many sub epitopes needed to be matched to a lot of proteins, so to do this within reasonable time, it utilizes an Aho-Corasick algorithm, this is a bit complicated but it sorta created a tree-structure where each branch is a sub epitope and if two sub-epitopes share a prefix they will share a branch until the differ where they branch off. there is more to it if you want to know more look at internet. Then the matching takes place where all sub-epitopes is compared to all protein, creating a lot of possible redundancy. We only want the longest match, meaning we only want 1 match per epitope, since one epitope might get matches on multiple of it sub epitopes, so we only keep the longest match per epitope-pathogen-protein pair. We also keep the pathogen proteins that did not match but flag them.
* Input: output from 3\_Protein\_meta\_data (pathogen protein table) and output from 1\_IEDB\_Wrangling.py (Autoimmune epitope table)
* Output: Table where the rows is of each pathogen-protein-epitope match or pathogen-protein unmatched, with all metadata from the initial input tables.
