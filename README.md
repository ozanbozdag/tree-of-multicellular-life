# tree-of-multicellular-life
Overview

This R script (generate_tree.R) automates the creation of a phylogenetic cladogram (Figure 1) depicting transitions to multicellularity across the tree of life, including archaeal, bacterial, and eukaryotic lineages. It was developed for the manuscript "Overcoming key challenges in the evolutionary transition to multicellularity" (under review).

Purpose

Use case: Generate a high-resolution cladogram for Figure 1, showing multicellular lineages and their unicellular relatives.

Data sources:

Taxonomic validation via NCBI Taxonomy Browser

Tree topology from the Open Tree of Life API

Dependencies

R packages: httr, jsonlite, rotl, ape, phytools

Internet connection to query the Open Tree of Life server

Script Workflow

Define species list: Specify target multicellular taxa (and close unicellular relatives) spanning Archaea, Bacteria, and Eukarya.

Resolve OTT IDs: Match species names to Open Tree Taxonomy (OTT) identifiers using tnrs_match_names().

Filter matches: Retain taxa with valid OTT IDs and sufficient match confidence.

Retrieve subtree: Fetch the induced phylogenetic subtree from the Open Tree of Life via tol_induced_subtree().

Polytomy resolution: Convert the tree to a fully bifurcating (binary) topology with multi2di().

Visualization & export:

Generate fan‐shaped and rectangular tree plots (PDF, TIFF, PNG) with customized fonts and margins.

Save Newick and Nexus tree files.

Summary output: Print file locations, tree statistics (number of tips/nodes), and sorted tip labels in the console.

Script Usage

# From your R console or script runner:
source("generate_tree.R")

Outputs are written to the tree_output/ directory under your working directory.

Output Files

multicellular_lineages_fan.pdf / .tiff / .png — fan‐layout cladogram

multicellular_lineages_rectangular.pdf — rectangular‐layout cladogram

multicellular_lineages.nwk — Newick format tree

multicellular_lineages.nex — Nexus format tree

Author: G. Ozan BozdagScripts & data sources documented within.
