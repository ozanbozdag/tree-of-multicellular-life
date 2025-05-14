#Author: G Ozan Bozdag
#Taxanomic validation:https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi
#Evolutionary tree source data:https://tree.opentreeoflife.org 

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
  
  # Defulfobacteria Clade
  "Magnetoglobus multicellularis",
  "Algorimarina butyrica",
  
  # Thermodesulfobacteriota Clade
  "Electronema",
  "Desulfovibrio vulgaris", # basionym: Nitratidesulfovibrio vulgaris
  
  #  Clostridia-Eubacteriales Clade
  "Candidatus Arthromitus",
  "Heliophilum fasciatum",
  
  # Bacillati-chloroflexota Clade
  "Chloroflexus",
  "Dehalococcoides mccartyi",
  
  # Pseudomonadati-proteobacteria Clade
  "Thioploca araucae",
  "Thiobacillus denitrificans",
  
  # Pseudomonadati-Myxococcota Clade
  "Anaeromyxobacter dehalogenans",
  "Chondromyces crocatus",
  
  # Chlorophyta+Streptophyta+Rhodophyta (confirmed that the output of this part of 
  # the tree topology matches the manually generated tree
  # in Umen&Herron 2021 Ann Rev Fig1):
  
  # 1  Chlorophyta-Ulvales Clade (Class: Ulvophyceae)
  "Ulva lactuca",
  "Oltmannsiellopsis viridis",
  
  # 2 Volvocine algae Clade
  "Chlamydomonas reinhardtii",
  # WARNING: The data on open tree of life server misplaces 
  # "Tetrabaena socialis" (appears as an outgroup which is incorrect), 
  # to manually correct for that, used Gonium (which is correctly placed within the clade), 
  # and rename it as "Tetrabaena socialis" on the PDF file:
  "Gonium multicoccum", #rename as "Tetrabaena socialis" on the PDF
  "Volvox carteri",
  
  # 3  Prasiolales Clade
  "Prasiola calophylla",
  
  # 4 Streptophyta-Charophytes Clade
  "Chara braunii", # Charophytes
  "Mesostigma viride", #  Streptophyta
  "Rosa gallica", #Land plants
  "Mesotaenium endlicherianum", # Zygnematophyceae
  
  # 5 Rhodophyta / red algae:
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