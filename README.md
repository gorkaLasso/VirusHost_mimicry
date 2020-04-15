# Structural relationships between viruses and hosts

## Overview

## Table of content

### virus_viralProt_mapping.txt
Tab delimited text file describing the viral dataset, including the virus and corresponding viral proteins and hosts. This dataset was originaly obtained from the virushostDB (https://www.genome.jp/virushostdb/). Each line corresponds to a viral protein.

Columns:

Column 0	Gene/Protein description	Name for the corresponding viral gene and/or protein
Column 1	NCBI Identifier	Genomic and/or protein NCBI identifier
Column 2	Ids for other databases	Identifier for additional databases (e.g. GI, uniprotKb, Interpro)
Column 3	Internal VirProt id A	Unique internal identifier A for viral protein
Column 4	Internal VirProt id B	Internal identifier B for viral protein (two proteins 100% identical will share the same identifier)											
Column 5	Virus taxid	Virus taxid											
Column 6	Virus Name	Virus name											
Column 7	Virus family	Virus family											
Column 8	Nucleic acid type	Nucleic acid type for virus											
Column 9	Host TaxGroup Level 1	Broadest taxonomic division for hosts (Bacteria, Invertebrate and Bacteria, Invertebrate and Plant Fungi, Invertebrate & Vertebrate, Plant Fungi, Vertebrate)
Column 10	Host TaxGroup Level 2	The Vertebrate taxonomic division is broken down into mammals/non-mammals											
Column 11	Host TaxGroup Level 3	The taxonomic divisions containing mammals are further broken down into human/non-human											
Column 12Host taxid	Taxonomic identifier for  virus hosts		
Column 13Host name	Host names for virus											

### Structural relationships