#Author: G Ozan Bozdag
#Taxanomic validation: https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi
#Evolutionary tree source data: https://tree.opentreeoflife.org 
#
# REFERENCE NOTES FOR MANUAL TOPOLOGY CORRECTIONS
#
# The Open Tree of Life (OTL) induced subtree contains several topological 
# inaccuracies that require manual correction in the final PDF. Corrections 
# are based on the following published phylogenies:
#
#  - BACTERIA: Hug et al. 2016, "A new view of the tree of life", 
#    Nature Microbiology, DOI: 10.1038/NMICROBIOL.2016.48
#    Used for: Cyanobacteria placement, Chloroflexi/Firmicutes grouping, 
#    separation of former "Deltaproteobacteria" (Myxococcota + Desulfobacterota) 
#    from Proteobacteria s.s. (Alpha/Beta/Gamma). Note: Hug et al. show that 
#    Proteobacteria is not monophyletic — Delta branches away from the others.
#
#  - HOLOZOA (Choanoflagellata, Filasterea, Ichthyosporea, Metazoa): 
#    Ruiz-Trillo et al. 2023 used for overall holozoan topology.
#    Brunet et al. 2019 used for choanoflagellate tip ordering 
#    (Salpingoeca rosetta as outgroup to Monosiga + Choanoeca).
#    Note: Choanoeca flexa not in OTL DB — used Choanoeca perplexa as 
#    placeholder, rename in final PDF. Chromosphaera perkinsii not in OTL DB 
#    — used Pirum gemmata as placeholder, rename in final PDF.
#
#  - GREEN ALGAE (Chlorophyta + Streptophyta + Rhodophyta): 
#    Umen & Herron 2021 (Annu. Rev. Genet.) Fig. 1 used to confirm topology. 
#    Note: OTL misplaces Tetrabaena socialis as an outgroup; used Gonium 
#    multicoccum (correctly placed) as placeholder, rename in final PDF.
# =============================================================================

# Load required libraries
library(httr)
library(jsonlite)
library(rotl)
library(ape)
library(phytools)

# Print working directory at start
cat("\nCurrent working directory:", getwd(), "\n")

# Define multicellular species list
species <- c(
  # Methanobacteriota-Archaea Clade
  "Methanosarcina barkeri", #couldn't find a unicellular relative distant 
  #from Halobacterial clade below
  
  #  Halobacteria-Archaea Clade
  "Haloferax volcanii",
  "Haloarcula marismortui",
  
  # Desulfobacterota Clade (former "Deltaproteobacteria" per Hug et al. 2016)
  # NOTE: manually group Magnetoglobus+Algorimarina and Electronema+Desulfovibrio
  # with Myxococcota (below) on a deep former-"Deltaproteobacteria" branch,
  # SEPARATE from Proteobacteria s.s. (Thioploca+Thiobacillus).
  "Magnetoglobus multicellularis",
  "Algorimarina butyrica",
  
  # Thermodesulfobacteriota Clade
  "Electronema",
  "Desulfovibrio vulgaris", # basionym: Nitratidesulfovibrio vulgaris
  
  #  Clostridia-Eubacteriales Clade (Firmicutes)
  # NOTE: manually group Arthromitus+Heliophilum as sister pair (Firmicutes),
  # then sister to Chloroflexi (below) per Hug et al. 2016.
  # OTL incorrectly nests Heliophilum/Arthromitus with Dehalococcoides.
  "Candidatus Arthromitus",
  "Heliophilum fasciatum",
  
  # Bacillati-chloroflexota Clade (Chloroflexi)
  # NOTE: manually group Chloroflexus+Dehalococcoides as sister pair,
  # then sister to Firmicutes (above) per Hug et al. 2016.
  "Chloroflexus",
  "Dehalococcoides mccartyi",
  
  # Pseudomonadati-proteobacteria Clade (Proteobacteria s.s.)
  # NOTE: keep Thioploca (Gamma) + Thiobacillus (Beta) together as
  # Proteobacteria s.s. clade, SEPARATE from former "Deltaproteobacteria"
  # (Myxococcota + Desulfobacterota above) per Hug et al. 2016 — these are
  # NOT sister groups; Delta branches away from Proteobacteria.
  "Thioploca araucae",
  "Thiobacillus denitrificans",
  
  # Pseudomonadati-Myxococcota Clade (former "Deltaproteobacteria" per Hug et al.)
  # NOTE: manually group Anaeromyxobacter+Chondromyces with Desulfobacterota
  # (above) on the deep former-"Deltaproteobacteria" branch; pull OUT from
  # any nesting with Proteobacteria s.s.
  "Anaeromyxobacter dehalogenans",
  "Chondromyces crocatus",
  
  # Cyanobacteria Clade
  # NOTE: OTL incorrectly nests this pair with Chloroflexi. Manually
  # reposition Anabaena+Gloeobacter as its own branch BEFORE (sister to)
  # the Chloroflexi+Firmicutes cluster, per Hug et al. 2016 rooted tree.
  "Anabaena sphaerica", # multicellular filamentous cyanobacterium with heterocysts
  "Gloeobacter violaceus", # unicellular, deepest-branching cyanobacterium (lacks thylakoids)
  
  # Chlorophyta+Streptophyta+Rhodophyta (confirmed that the output of this part of 
  # the tree topology matches the manually generated tree
  # in Umen&Herron 2021 Ann Rev Fig1):
  
  # Chlorophyta-Ulvales Clade (Class: Ulvophyceae)
  "Ulva lactuca",
  "Oltmannsiellopsis viridis",
  
  # Volvocine algae Clade
  "Chlamydomonas reinhardtii",
  # WARNING: The data on open tree of life server misplaces 
  # "Tetrabaena socialis" (appears as an outgroup which is incorrect), 
  # to manually correct for that, used Gonium (which is correctly placed within the clade), 
  # and rename it as "Tetrabaena socialis" on the PDF file:
  "Gonium multicoccum", #rename as "Tetrabaena socialis" on the PDF
  "Volvox carteri",
  
  # Prasiolales Clade
  "Prasiola calophylla",
  
  # Streptophyta-Charophytes Clade
  "Chara braunii", # Charophytes
  "Mesostigma viride", #  Streptophyta
  "Rosa gallica", #Land plants
  "Mesotaenium endlicherianum", # Zygnematophyceae
  
  # hodophyta / red algae:
  "Pyropia tenera",
  "Porphyridium purpureum",
  
  #Holozoa (needs to be manually 
  #updated, used the tree in Ruiz-Trillo et al. 2023):
  "Pan paniscus", #Metazoa
  
  #Choano (note: need to switch perplexa/flexa and rosetta's positions
  # at the tip of the tree, following Brunet et al 2019):
  "Salpingoeca rosetta", #make this outgroup to the two below after generating the tree
  "Monosiga brevicollis",
  "Choanoeca perplexa", # Choanoeca flexa, not found in the OTL database,
  #rename Perplexa as Choanoeca flexa after generating the tree
  
  #Filasterea:
  "Capsaspora owczarzaki",
  "Ministeria vibrans",
  
  #Ichthyosporea:
  "Sphaeroforma arctica", 
  "Pirum gemmata", # C. perkinsii, not found in the OTL database,
  #rename Pirum as Chromosphaera perkinsii after generating the tree
  
  # Lentinula
  "Lentinula edodes",
  "Rhodotorula mucilaginosa",
  
  # Aspergillus
  "Aspergillus fumigatus",
  "Saccharomyces cerevisiae",
  
  # Nereocystis luetkeana Clade
  "Nereocystis luetkeana",
  "Chattonella subsalsa",
  
  # Synura uvella Clade
  "Synura uvella",
  "Mallomonas akrokomos",
  
  # Albugo candida Clade
  "Albugo candida",
  "Hyphochytrium catenoides",
  
  # Viridiuvalis adhaerens Clade
  "Viridiuvalis adhaerens",
  "Lotharella globosa",
  
  # Zoothamnium niveum Clade
  "Zoothamnium niveum",
  "Epistylis anastatica",
  
  # Dictyostelium discoideum Clade
  "Dictyostelium discoideum",
  "Dermamoeba algensis",
  
  # Copromyxa protea Clade
  "Copromyxa protea",
  "Amoeba proteus",
  
  # Fonticula alba Clade
  "Fonticula alba",
  "Nuclearia simplex", 
  
  # Sorogena stoianovitchae Clade
  "Sorogena stoianovitchae",
  "Cyrtolophosidida",
  
  # Sorodiplophrys stercorea Clade
  "Sorodiplophrys stercorea",
  "Stellarchytrium dubum",
  
  # Guttulinopsis nivea Clade
  "Guttulinopsis nivea",
  
  # Acrasis rosea Clade
  "Acrasis rosea",
  "Naegleria gruberi"
)

# Match names to OTT IDs
cat("\nAttempting to match species names to OTT IDs...\n")
taxa <- rotl::tnrs_match_names(species)
print("Match results:")
print(taxa)

# Print specific information about Arthromitus
arthromitus_match <- taxa[taxa$search_string == "Arthromitus", ]
cat("\nSpecific match information for Arthromitus:\n")
print(arthromitus_match)

# Filter for valid taxa with good matches
valid_taxa <- taxa[!is.na(taxa$ott_id) & taxa$score >= 0.5, ]  # Lowered threshold to 0.5 to include more matches

# Print any excluded taxa
excluded <- taxa[!taxa$search_string %in% valid_taxa$search_string, ]
if(nrow(excluded) > 0) {
  cat("\nWarning: The following taxa were excluded:\n")
  print(excluded)
}

# Generate tree
cat("\nGenerating phylogenetic tree...\n")
tr <- rotl::tol_induced_subtree(ott_id(valid_taxa), label="name")

# Convert polytomies to binary trees (nodes with exactly two descendants)
cat("\nConverting polytomies to dichotomies...\n")
tr_binary <- ape::multi2di(tr, random=TRUE)

# Create output directory if it doesn't exist
output_dir <- file.path(getwd(), "tree_output")
dir.create(output_dir, showWarnings = FALSE)

# Save files with full paths
fan_pdf_file <- file.path(output_dir, "multicellular_lineages_fan.pdf")
rect_pdf_file <- file.path(output_dir, "multicellular_lineages_rectangular.pdf")
fan_tiff_file <- file.path(output_dir, "multicellular_lineages_fan.tiff")
fan_png_file <- file.path(output_dir, "multicellular_lineages_fan.png")
newick_file <- file.path(output_dir, "multicellular_lineages.nwk")
nexus_file <- file.path(output_dir, "multicellular_lineages.nex")

# IMPROVED FAN LAYOUT with SMALLER FONT - FIXED FOR PDF
# Key parameter: xpd=TRUE allows plotting outside figure region
pdf(fan_pdf_file, width=20, height=20, useDingbats=FALSE)
par(mar=c(3,3,3,3), xpd=TRUE)  # xpd=TRUE allows plotting outside the figure region

# Plot with smaller font size
plot.phylo(tr_binary, 
           type="fan",           
           cex=0.65,            # REDUCED font size
           font=1,              # Regular font (not bold)
           label.offset=0.05,   # Keep labels close to tips
           show.tip.label=TRUE, # Force display of tip labels
           show.node.label=FALSE,
           edge.width=1.0,
           tip.color="black",   
           align.tip.label=FALSE)

# Title and legend with smaller font
title("Phylogenetic Tree of Multicellular Lineages (Fan Layout)\nBased on Open Tree of Life", 
      cex.main=1.2)      # Smaller title
legend("bottomleft", 
       legend=paste("Number of tips:", length(tr_binary$tip.label)),
       bty="n",
       cex=0.9)          # Smaller legend
dev.off()

# High-resolution TIFF with smaller font - also with xpd=TRUE
tiff(fan_tiff_file, width=2400, height=2400, res=300)
par(mar=c(3,3,3,3), xpd=TRUE)  # xpd=TRUE allows plotting outside the figure region

# Plot with smaller font for TIFF
plot.phylo(tr_binary, 
           type="fan",           
           cex=0.7,              # Slightly larger than PDF but still small
           font=1,               # Regular font
           label.offset=0.05,   
           show.tip.label=TRUE,  
           show.node.label=FALSE,
           edge.width=1.2,
           tip.color="black",    
           align.tip.label=FALSE) 

title("Phylogenetic Tree of Multicellular Lineages (Fan Layout)\nBased on Open Tree of Life", 
      cex.main=1.2)
legend("bottomleft", 
       legend=paste("Number of tips:", length(tr_binary$tip.label)),
       bty="n",
       cex=0.9)
dev.off()

# Rectangular layout with xpd=TRUE
pdf(rect_pdf_file, width=25, height=18)
par(mar=c(3,12,3,3), xpd=TRUE)  # xpd=TRUE allows plotting outside the figure region

plot.phylo(tr_binary,
           type="phylogram",
           direction="rightwards",
           cex=0.8,                # Moderate font size
           font=1,                 # Regular font
           label.offset=0.1,        
           show.tip.label=TRUE,     
           show.node.label=FALSE,
           edge.width=1.2,
           tip.color="black",      
           align.tip.label=TRUE)    

title("Phylogenetic Tree of Multicellular Lineages (Rectangular Layout)", 
      cex.main=1.2,
      line=0)
dev.off()

# Save other formats
write.tree(tr, file=newick_file)
write.nexus(tr, file=nexus_file)

# Print file locations and summary
cat("\nFiles have been saved to:\n")
cat("Fan PDF file:", fan_pdf_file, "\n")
cat("Rectangular PDF file:", rect_pdf_file, "\n")
cat("Fan TIFF file (high resolution):", fan_tiff_file, "\n")
cat("Newick file:", newick_file, "\n")
cat("NEXUS file:", nexus_file, "\n")

# Plot in RStudio window with improved visibility
par(mar=c(3,3,3,3), xpd=TRUE)  # xpd=TRUE allows plotting outside the figure region
plot.phylo(tr_binary, 
           type="fan",           
           cex=0.65,             # Same small font size as PDF
           font=1,              
           label.offset=0.05,   
           show.tip.label=TRUE, 
           show.node.label=FALSE,
           edge.width=1.0,
           tip.color="black",   
           align.tip.label=FALSE)
title("Phylogenetic Tree of Multicellular Lineages\nBased on Open Tree of Life", 
      cex.main=1.0)

# SAVE THE CURRENT PLOT DIRECTLY FROM RSTUDIO
cat("\nTo save the current plot from RStudio, run these lines:\n")
cat('dev.copy(png, "', fan_png_file, '", width=1500, height=1500, res=150)', "\n", sep="")
cat("dev.off()", "\n")

# Print tree statistics
cat("\nTree Statistics:\n")
cat("Number of tips:", length(tr$tip.label), "\n")
cat("Number of internal nodes:", tr$Nnode, "\n")
cat("\nTip labels:\n")
print(sort(tr$tip.label))