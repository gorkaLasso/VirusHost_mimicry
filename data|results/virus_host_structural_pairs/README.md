# Virus-Host structural pairs

Tab delimited text files describing structural relationships between viral proteins and proteins in non-viral organisms. Each file within the directory corresponds to a set of viruses grouped together according to the taxonomic division of their host.

Important! The run for SARS-CoV-2 proteins was carried out using an updated version of the PDB. Therefore, human proteins mimicking SARS-CoV-2 but not other coronaviruses might be the result of having used a more up-to-date database.

     Column 0   Internal_virprotID_A: Unique internal identifier A for viral protein

     Column 1   str_temp: Intermediate structural template used to infer a structural relationship between proteins

     Column 2   strl_neig_uni: Uniprot accesion code for structural neighbor

     Column 3   str_neigh_pdb: PDB code for structural neighbor

     Column 4   sas: Structural Alignment Score

     Column 5   seqHomology: Significant sequence homology between viral and structural neighbor (seq homolog if e-value < 1x10^-6)
