# Data & Results

## Table of contents
* [Mapping](#mapping)
* [Human_infecting_viral_families](#human_infecting_viral_families)
* [Structural_relationships](#structural_relationships)


## Mapping
* mapping/virus_viralProt_mapping.txt

Tab delimited text file describing the viral dataset, including the virus and corresponding viral proteins and hosts. This dataset was originaly obtained from the virushostDB (https://www.genome.jp/virushostdb/). Each line corresponds to a viral protein. 

     Column 0   Gene/Protein description: Name for the corresponding viral gene and/or protein

     Column 1   NCBI Identifier: Genomic and/or protein NCBI identifier

     Column 2   Ids for other databases: Identifier for additional databases (e.g. GI, uniprotKb, Interpro)

     Column 3   Internal VirProt id A: Unique internal identifier A for viral protein

     Column 4   Internal VirProt id B: Internal identifier B for viral protein (two proteins 100% identical will share the same identifier)

     Column 5   Virus taxid: virus taxonomic identifier

     Column 6   Virus Name

     Column 7   Virus family

     Column 8   Nucleic acid type

     Column 9   Host TaxGroup Level 1: Broadest taxonomic division for hosts (Bacteria, Invertebrate and Bacteria, Invertebrate and Plant Fungi, Invertebrate & Vertebrate, Plant Fungi, Vertebrate)

     Column 10  Host TaxGroup Level 2: The Vertebrate taxonomic division is broken down into mammals/non-mammals

     Column 11  Host TaxGroup Level 3: The taxonomic divisions containing mammals are further broken down into human/non-human

     Column 12  Host taxid: Taxonomic identifier for virus hosts

     Column 13  Host name: Names for virus hosts

## Human_infecting_viral_families
* human_infecting_viral_families/*.txt file

Tab delimited text files (one file per viral family) describing the percentage of human-infecting viruses within a viral family that structurally mimic a particular human protein (least one viral protein structurally mimics the human protein)

     Column 0   Human protein (Uniprot AC) mimicking >= viruses within the viral family

     Column 1   Number of viruses within the family structurally mimicking the human protein

     Column 2   Total number of viruses within a viral family

     Column 3   Percentages of viruses within the family structurally mimicking the human protein

## Structural_relationships
* virus_host_structural_pairs/*.txt files

Tab delimited text files describing structural relationships between viral proteins and proteins in non-viral organisms. Each file within the directory corresponds to a set of viruses grouped together according to the taxonomic division of their host.

Important! The run for SARS-CoV-2 proteins was carried out using an updated version of the PDB. Therefore, human proteins mimicking SARS-CoV-2 but not other coronaviruses might be the result of having used a more up-to-date database.

     Column 0   Internal_virprotID_A: Unique internal identifier A for viral protein

     Column 1   str_temp: Intermediate structural template used to infer a structural relationship between proteins

     Column 2   strl_neig_uni: Uniprot accesion code for structural neighbor

     Column 3   str_neigh_pdb: PDB code for structural neighbor

     Column 4   sas: Structural Alignment Score

     Column 5   seqHomology: Significant sequence homology between viral and structural neighbor (seq homolog if e-value < 1x10^-6)