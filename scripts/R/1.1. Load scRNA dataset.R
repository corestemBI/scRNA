# Corestem project - single_cell_practice.proj
# 20220125 Sehwan Chun at Corestem, Inc.
# 1.1. Load scRNA dataset.R

#Before Running scripts 
options(java.parameters = c("-XX:+UseConcMarkSweepGC", "-Xmx8192m"))
memory.limit(size = 100000000000)
source("./functions.R")

#### 1. Data loading ####

# Load the MSC data sets
bm_msc  = Read10X(data.dir = "X:/Single_cell/HN00160068_10X_CellRanger_Outs/run_count_BM-MSC/outs/filtered_feature_bc_matrix/")
ucb_msc = Read10X(data.dir = "X:/Single_cell/HN00160068_10X_CellRanger_Outs/run_count_UCB-MSC/outs/filtered_feature_bc_matrix/")
wj_msc  = Read10X(data.dir = "X:/Single_cell/HN00160068_10X_CellRanger_Outs/run_count_WJ-MSC/outs/filtered_feature_bc_matrix/")
cpj_msc = Read10X(data.dir = "X:/Single_cell/HN00160068_10X_CellRanger_Outs/run_count_CPJ-MSC/outs/filtered_feature_bc_matrix/")

# Make seurat object from 10X data
bm_msc  = run_seurat_object(bm_msc,  "BM-MSC")
ucb_msc = run_seurat_object(ucb_msc, "UCB-MSC")
wj_msc  = run_seurat_object(wj_msc,  "WJ-MSC")
cpj_msc = run_seurat_object(cpj_msc, "CPJ-MSC")

#### 2. Data integrating by SCT method ####

# Get common features of MSCs and integration through SCT
all_msc_list = c(bm_msc, ucb_msc, wj_msc, cpj_msc)
msc_all_combined = run_seurat_integrated(all_msc_list)
rm(bm_msc);rm(cpj_msc);rm(ucb_msc);rm(wj_msc);rm(all_msc_list) # For memory

#### 3. save and load ####

# Save as h5seurat image
SaveH5Seurat(msc_all_combined, "X:/single_cell/HN00160068_10X_Rimages_Outs/seurat_combined_normal.image.h5seurat", overwrite = T, verbose = T)
msc_all_combined = LoadH5Seurat("X:/single_cell/HN00160068_10X_Rimages_Outs/seurat_combined_normal.image.h5seurat")
