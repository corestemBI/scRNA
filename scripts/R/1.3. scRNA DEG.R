# Corestem project - single_cell_practice.proj
# 20220125 Sehwan Chun at Corestem, Inc.
# 1.3. scRNA DEG.R

#Before Running scripts 
options(java.parameters = c("-XX:+UseConcMarkSweepGC", "-Xmx8192m"))
memory.limit(size = 100000000000)
source("./functions.R")

#### 1. Load Seurat Object ####

msc_all_combined = LoadH5Seurat("X:/single_cell/HN00160068_10X_Rimages_Outs/seurat_combined_normal.image.h5seurat")

#### 2. Subtype-specific markers ####

# Get MSC subtype-specific feature genes
msc_all_combined = run_seurat_normalization(msc_all_combined)
msc_all_markerslist = run_seurat_markers(objects = msc_all_combined, project_name = "MSC_ALL")
msc_all_markerslist2 = run_seurat_markers_all(objects = msc_all_combined, project_name = "MSC_ALL")
rm(list = grep("tmp", ls(),value = T))
run_seurat_ridge_markers(msc_all_combined, "MSC_ALL")
run_seurat_violin_markers(msc_all_combined, "MSC_ALL")

# DEG to pathway 
run_DEG_profiler(DEGs = "./MSC_ALL.xlsx", sheets = 1, project_name = "MSC_ALL_res01")
run_DEG_profiler(DEGs = "./MSC_ALL.xlsx", sheets = 2, project_name = "MSC_ALL_res06")
run_DEG_profiler(DEGs = "./MSC_ALL.xlsx", sheets = 3, project_name = "MSC_ALL_origin")
run_DEG_multiprofiler(DEGs = "./MSC_ALL.xlsx", sheets = 1, project_name = "MSC_ALL_res01")
run_DEG_multiprofiler(DEGs = "./MSC_ALL.xlsx", sheets = 2, project_name = "MSC_ALL_res06")
run_DEG_multiprofiler(DEGs = "./MSC_ALL.xlsx", sheets = 3, project_name = "MSC_ALL_origin")






