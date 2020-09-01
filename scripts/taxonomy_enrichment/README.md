# taxonomy_enrichment

Taxonomic enrichment for a group of viruses using the set of proteins that are structurally similar to the viral proteins. Viruses are grouped according to the taxonomic division of their host. Taxonomic division enrichment was computed with a hypergeometric test that describes the significance of having k structural neighbors belonging to a particular taxonomic division (out of n total structural neighbors for a group of viruses) given the entire set of structurally solved proteins in the PDB of size N that contains K proteins from the same taxonomic division. To minimize experimental bias of multiple PDB entries for the same protein, structurally solved proteins were mapped to their Uniprot accession codes using SIFTS mapping (https://www.ebi.ac.uk/pdbe/docs/sifts/).

## Table of contents
* [Directories](#Directories)
* [1_identify_templates.pl](#1_identify_templates.pl)
* [2_identify_neighbors.pl](#2_identify_neighbors.pl)
* [3_hypergeometric.pl](#3_hypergeometric.pl)
* [4_identity_seqBasedTaxClusters.pl](#4_identity_seqBasedTaxClusters.pl)
* [5_hypergeometric_seqBasedTaxClusters.pl](#5_hypergeometric_seqBasedTaxClusters.pl)
* [enrichment_sum_4_manuscript.xlsx](#enrichment_sum_4_manuscript.xlsx)


## Directories (Bacteria, Human, Invertebrate, PlantFungi, Vertebrate)

Taxonomy enrichment results for Bacteria-infecting viruses, Human-infecting viruses, Invertebrate-infecting viruses, PlantFungi-infecting viruses and Vertebrate-infecting viruses. Files contained within each directory:

* enrichment_2.5.txt

     Hypergeometric test using all structural mimics

* enrichment_subCluster40_2.5.txt

     Hypergeometric test using a set of structural mimics filtered at 40% sequence identity

* representative_template_space.txt

     Set of structural templates used to infer structurally similar proteins (pdb, pdb_start, pdb_end, uniprot AC, taxonomic division)

* subCluster40_neighbor_space_2.5.txt

     Set of inferred structural mimics filtered at 40% sequence identity

* uniprot_neighbor_space.txt

     Set of inferred structural mimics (uniprot, host taxonomic ID, host name, host taxonomic division, lowest structural alignment score)

## 1_identify_templates.pl

## 2_identify_neighbors.pl

## 3_hypergeometric.pl

## 3_hypergeometric_taxID.pl

## 4_identity_seqBasedTaxClusters.pl

## 5_hypergeometric_seqBasedTaxClusters.pl

## enrichment_sum_4_manuscript.xlsx
