# Corestem project - single_cell_practice.proj
# 20220125 Sehwan Chun at Corestem, Inc.
# 1.2. scRNA characteristics and clustering.R

#Before Running scripts 
options(java.parameters = c("-XX:+UseConcMarkSweepGC", "-Xmx8192m"))
memory.limit(size = 100000000000)
source("./functions.R")

#### 1. Load Seurat Object ####

msc_all_combined = LoadH5Seurat("X:/single_cell/HN00160068_10X_Rimages_Outs/seurat_combined_normal.image.h5seurat")

#### 2. Expression strength of lineage-related genes and effector molecules ####

# Expression strength of lineage-related genes
msc_all_combined = run_seurat_lineage_score(msc_all_combined)
msc_all_combined = run_seurat_lineage_assign(msc_all_combined)

# Expression strength of ALS/MSC effector molecules
msc_all_combined = run_seurat_ALSMSC_score(msc_all_combined)

#### 3. Dimension reduction ####

# Dimplot using UMAP & t-SNE algorithm (Origin)
msc_all_combined = run_seurat_umap(msc_all_combined, "MSC_ALL_01", 0.1)
msc_all_combined = run_seurat_tsne(msc_all_combined, "MSC_ALL_01", 0.1)
msc_all_combined = run_seurat_umap(msc_all_combined, "MSC_ALL_06", 0.6)
msc_all_combined = run_seurat_tsne(msc_all_combined, "MSC_ALL_06", 0.6)

msc_all_combined = run_seurat_umap_subset(msc_all_combined, "MSC_ALL_01", 0.1)
msc_all_combined = run_seurat_tsne_subset(msc_all_combined, "MSC_ALL_01", 0.1)
msc_all_combined = run_seurat_umap_subset(msc_all_combined, "MSC_ALL_06", 0.6)
msc_all_combined = run_seurat_tsne_subset(msc_all_combined, "MSC_ALL_06", 0.6)

# Dimplot using UMAP & t-SNE algorithm (Lineage)
run_seurat_umap_lineage(msc_all_combined, "MSC_lineage_01", 0.1)
run_seurat_tsne_lineage(msc_all_combined, "MSC_lineage_01", 0.1)

# Dimplot using UMAP & t-SNE algorithm (score from ALS/MSC effector gene module)
run_seurat_umap_effector_score(msc_all_combined, "MSC_effector_01", 0.1)
run_seurat_tsne_effector_score(msc_all_combined, "MSC_effector_01", 0.1)

# imbalance of MSC distribution with chisq test
run_seurat_lineage_chisq(msc_all_combined, "MSC_ALL") # Lineage-Origin
run_seurat_cluster_source_chisq(msc_all_combined, "MSC_ALL")
run_seurat_cluster_lineage_chondrogenic_chisq(msc_all_combined, "MSC_ALL")
run_seurat_cluster_lineage_adipogenic_chisq(msc_all_combined, "MSC_ALL")
run_seurat_cluster_lineage_osteogenic_chisq(msc_all_combined, "MSC_ALL")

# Boxplot (effector scores/Lineage ~ Origin/Clusters)
run_seurat_msc_score_boxplot(msc_all_combined, "MSC_effector_01")





