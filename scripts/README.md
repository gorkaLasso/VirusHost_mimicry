# Scripts

## Table of contents
* [pathway_enrichment_humanVirFam](#pathway_enrichment_humanVirFam)
* [taxonomy_enrichment](#taxonomy_enrichment)

## pathway_enrichment_humanVirFam

For each viral family, the space of structural neighbors was used for enrichment analysis. Enrichment of biological pathways was determined using David, background corrected using the 5,841 human protein structures extracted from PDB. A Bonferroni corrected p-value < 0.01 was used to identify enriched biological ontologies. Viral families and biological pathways enriched in at least one viral family were clustered using R and the heatmap.2 function within gplots package (values in original matrix: -log10 P-valueBonferroni; distance metric: Euclidean, method: complete). 

## taxonomy_enrichment

Taxonomic enrichment for a group of viruses using the set of proteins that are structurally similar to the viral proteins. Viruses are grouped according to the taxonomic division of their host. Taxonomic division enrichment was computed with a hypergeometric test that describes the significance of having k structural neighbors belonging to a particular taxonomic division (out of n total structural neighbors for a group of viruses) given the entire set of structurally solved proteins in the PDB of size N that contains K proteins from the same taxonomic division. To minimize experimental bias of multiple PDB entries for the same protein, structurally solved proteins were mapped to their Uniprot accession codes using SIFTS mapping (https://www.ebi.ac.uk/pdbe/docs/sifts/).

