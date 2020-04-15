# Structural relationships between viruses and hosts

## Table of contents
* [Overview](#overview)
* [Viral_dataset](#viral_dataset)
* [Structural_relationships](#structural_relationships)

## Overview
We employ sequence-based methods to identify proteins that have similar structures to queried viral proteins and then use structural alignment to find “structural neighbors” of viral proteins. We applied the approach to a set of 337,493 viral proteins representing 7,486 viruses across a broad host taxonomic range, including bacteria, plants and fungi, invertebrates and vertebrates. Our survey identified over 6,000,000 structural relationships between proteins in viruses and non-viral organisms.

## Viral_dataset
* virus_viralProt_mapping.txt

Tab delimited text file describing the viral dataset, including the virus and corresponding viral proteins and hosts. This dataset was originaly obtained from the virushostDB (https://www.genome.jp/virushostdb/). Each line corresponds to a viral protein.

Columns:

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

## Structural relationships
* str_pairs/*.txt files

Tab delimited text file describing the inferred structural relationships inferred through an intermediate structural template.


     Column 0   Internal_virprotID_A: Unique internal identifier A for viral protein

     Column 1   str_temp: Intermediate structural template used to infer a structural relationship between proteins

     Column 2   strl_neig_uni: Uniprot accesion code for structural neighbor

     Column 3   str_neigh_pdb: PDB code for structural neighbor

     Column 4   sas: Structural Alignment Score

     Column 5   seqHomology: Significant sequence homology between viral and structural neighbor (seq homolog if e-value < 1x10^-6)
