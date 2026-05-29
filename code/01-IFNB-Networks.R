#~~~~~~~~~~~~~~~~~~~~~~~~ README ~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 
# Title: 01-IFNB-Networks.R
# Description: Demonstrate network analysis for pathway data
#
# Author: Mike Martinez
# Lab: Compute and Conquer
# Project: May 2026 Blog
# Date created: May 28th, 2026
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# LOAD LIBRARIES AND SET PATHS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#----- Libraries
library(Seurat)
library(SeuratData)
library(GDSCtools)
library(org.Hs.eg.db)
library(clusterProfiler)

#----- Set project directory
wd <- "/users/mike/Desktop/MayBlog/"

#----- Specify directories
outputDir <- paste0(wd, "outputs/")
figDir <- paste0(outputDir, "figures/")

#----- Create directories
dirs <- c(outputDir, figDir)
lapply(dirs, function(d) if (!dir.exists(d)) dir.create(d))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# LOAD IN THE IFNB DATASET
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#----- Install the dataset and load it
InstallData("ifnb")
seurat <- LoadData("ifnb")

#----- Set Default Assay to RNA
DefaultAssay(seurat) <- "RNA"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# DIFFERENTIAL EXPRESSION ANALYSIS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#----- Specify clusters
clusters <- unique(seurat$seurat_annotations)

#----- Run differential expresison analysis
de_list <- run_cluster_de(
  obj = seurat,
  clusters = clusters,
  cluster_col = "seurat_annotations",
  ident_col = "stim",
  ident1 = "STIM",
  ident2 = "CTRL")
write_de_results(de_list, paste0(outputDir, "Ctrl_v_Stim/DEGs/"), type = "DE", split_direction = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# OVER REPRESENTATION ANALYSIS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#----- Set universe
universe <- rownames(seurat)
set.seed(1234)

#----- GO-ORA analysis
go_list <- run_go_ora(de_list,
                      org_db = org.Hs.eg.db,
                      universe = universe,
                      outputDir = paste0(outputDir, "Ctrl_v_Stim/ORA/"))
write_de_results(go_list, paste0(outputDir, "Ctrl_v_Stim/ORA/"), type = "GO", split_direction = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# PATHWAY OVERLAP NETWORKS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#----- Network construction
networks <- run_comparison_networks(
  go_list    = go_list,
  de_list    = de_list,
  comp_label = "Ctrl_v_Stim",
  base_outdir = outputDir
)

#----- Extract the GO results with clustering
for (i in names(networks)) {
  x <- networks[[i]]$go_term_df
  cluster <- gsub("/", "_", i)
  name <- paste0(cluster, "_cluster_membership.csv")
  write.csv(x, file = paste0(outputDir, "Ctrl_v_Stim/networks/", name))
}





