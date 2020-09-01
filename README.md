# Structural relationships between viruses and hosts

## Table of contents
* [Overview](#Overview)
* [Data|results](#Data|Results)
* [Scripts](#Scripts)

## Overview
We employ sequence-based methods to identify proteins that have similar structures to queried viral proteins and then use structural alignment to find “structural neighbors” of viral proteins. We applied the approach to a set of 337,493 viral proteins representing 7,486 viruses across a broad host taxonomic range, including bacteria, plants and fungi, invertebrates and vertebrates. Our survey identified over 6,000,000 structural relationships between proteins in viruses and non-viral organisms.

## Data|results
* mapping 

    Mapping of viral proteins, viruses and their hosts as described by virushostDB (https://www.genome.jp/virushostdb/)

* virus_host_structural_pairs

    Inferred structural relationships between viral proteins and proteins in non-viral organisms

* human_infecting_viral_families

    Human structural mimics within human-infecting viral families

## Scripts
* taxonomy_enrichment

    Taxonomic enrichment for a group of viruses using the set of proteins that are structurally similar to the viral proteins. Viruses are grouped according to the taxonomic division of their host.
    
* pathway_enrichement_humanVirFam

    Clustering of the pathways enriched within a particular viral family (human-infecting viruses) based on the space of structural mimics.