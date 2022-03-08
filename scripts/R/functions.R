# Corestem project - single_cell_practice.proj
# 20220125 Sehwan Chun at Corestem, Inc.
# function.R

#### 1. Library Loading ####
packs = c("Seurat", "SeuratDisk", "dplyr", "patchwork", "ggpubr", "DoubletFinder", "xlsx", "ggstatsplot",
          "clusterProfiler", "org.Hs.eg.db", "AnnotationHub", "MAST", "ggsci", "RColorBrewer", "SeuratWrappers")
lapply(packs, require, character.only = TRUE)
rm(packs)

#### 2. Data loading and integrating ####
run_seurat_object = function(msc_data, project_name){
    
    msc_data  = CreateSeuratObject(counts = msc_data,  project = project_name, min.cells =3, min.features = 300)
    print(paste0("Raw Cells: ", ncol(msc_data)))
    # The [[ operator can add columns to object metadata. This is a great place to stash QC stats
    msc_data[["percent.mt"]]  =  PercentageFeatureSet(msc_data,  pattern = "^MT-")
    
    # Visualize QC metrics as a violin plot
    bqc_plot = VlnPlot(msc_data,  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
    
    # Remove doublet or low-quality QC cells! 
    msc_data = subset(msc_data,  subset = nFeature_RNA <= 6000 & nFeature_RNA >= 300 & percent.mt <= 5)
    print(paste0("Filtered Cells: ", ncol(msc_data)))
    aqc_plot = VlnPlot(msc_data,  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
    
    # Check cell cycle difference
    msc_data = SCTransform(msc_data,  do.scale = FALSE, do.center = FALSE, verbose = FALSE)
    msc_data = CellCycleScoring(msc_data,  s.features = cc.genes.updated.2019$s.genes,
                              g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = TRUE)
    
    ggsave(bqc_plot, filename = paste0("./outputs/",project_name,"_beforeQC.tiff"), width = 8, height = 8, dpi = 600,
           device = "tiff", compression = "lzw")
    ggsave(aqc_plot, filename = paste0("./outputs/",project_name,"_afterQC.tiff"),  width = 8, height = 8, dpi = 600,
           device = "tiff", compression = "lzw")
    ggsave(ggarrange(bqc_plot, aqc_plot, ncol = 2, labels = "auto"), width = 16, height = 8, dpi = 600,
           device = "tiff", compression = "lzw",
           filename = paste0("./outputs/",project_name,"_baQC.tiff"))
    
    # Get SCT count with vars.to.regress (MT, S, G2M)
    DefaultAssay(msc_data)  = "RNA"
    msc_data  = SCTransform(msc_data, vars.to.regress  = c("percent.mt", "S.Score", "G2M.Score"), verbose = FALSE)
    return(msc_data)
}
run_seurat_object_not_cc = function(msc_data, project_name){
    
    msc_data  = CreateSeuratObject(counts = msc_data,  project = project_name, min.cells =3, min.features = 300)
    print(paste0("Raw Cells: ", ncol(msc_data)))
    # The [[ operator can add columns to object metadata. This is a great place to stash QC stats
    msc_data[["percent.mt"]]  =  PercentageFeatureSet(msc_data,  pattern = "^MT-")
    
    # Visualize QC metrics as a violin plot
    bqc_plot = VlnPlot(msc_data,  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
    
    # Remove doublet or low-quality QC cells! 
    msc_data = subset(msc_data,  subset = nFeature_RNA <= 6000 & nFeature_RNA >= 300 & percent.mt <= 5)
    print(paste0("Filtered Cells: ", ncol(msc_data)))
    aqc_plot = VlnPlot(msc_data,  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
    
    # Check cell cycle difference
    msc_data = SCTransform(msc_data,  do.scale = FALSE, do.center = FALSE, verbose = FALSE)
    msc_data = CellCycleScoring(msc_data,  s.features = cc.genes.updated.2019$s.genes,
                                g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = TRUE)
    
    ggsave(bqc_plot, filename = paste0("./outputs/",project_name,"_beforeQC.tiff"), width = 8, height = 8, dpi = 600,
           device = "tiff", compression = "lzw")
    ggsave(aqc_plot, filename = paste0("./outputs/",project_name,"_afterQC.tiff"),  width = 8, height = 8, dpi = 600,
           device = "tiff", compression = "lzw")
    ggsave(ggarrange(bqc_plot, aqc_plot, ncol = 2, labels = "auto"), width = 16, height = 8, dpi = 600,
           device = "tiff", compression = "lzw",
           filename = paste0("./outputs/",project_name,"_baQC.tiff"))
    
    # Get SCT count with vars.to.regress (MT, S, G2M)
    DefaultAssay(msc_data)  = "RNA"
    msc_data  = SCTransform(msc_data, vars.to.regress  = c("percent.mt"), verbose = FALSE)
    return(msc_data)
}

run_seurat_integrated = function(object_list){
    msc_features = SelectIntegrationFeatures(object.list = object_list, nfeatures = 3000)
    integrated_msc = PrepSCTIntegration(object.list = object_list, anchor.features = msc_features)
    msc_anchors = FindIntegrationAnchors(object.list = integrated_msc, normalization.method = "SCT",
                                         anchor.features = msc_features)
    msc_combined_sct = IntegrateData(anchorset = msc_anchors, normalization.method = "SCT")
    DefaultAssay(msc_combined_sct) = "integrated"
    return(msc_combined_sct)
}

#### 3. Cell cycles ####
run_seurat_cellcyle = function(objects, project_name){
    
    Idents(objects) = objects@meta.data$Phase

    objects = RunPCA(objects, features = c(cc.genes.updated.2019$s.genes, cc.genes.updated.2019$g2m.genes))
    DimPlot(objects)
    
    p1 = DimPlot(objects)
    p1 = ggpar(p1, title = "")
    
    ggsave(p1, filename = paste0("./outputs/",project_name,"_cellcycle.tiff"), width = 8, height = 8, dpi = 600,
           device = "tiff", compression = "lzw")
} # Use before integration separately
run_seurat_cellcyle_subset = function(objects, project_name){
    
    Idents(objects) = objects@meta.data$orig.ident

    bm_cells  = WhichCells(objects, idents = "BM-MSC" )
    cpj_cells = WhichCells(objects, idents = "CPJ-MSC")
    ucb_cells = WhichCells(objects, idents = "UCB-MSC")
    wj_cells  = WhichCells(objects, idents = "WJ-MSC" )
    
    Idents(objects) = objects@meta.data$Phase
    
    objects = RunPCA(objects, features = c(cc.genes.updated.2019$s.genes, cc.genes.updated.2019$g2m.genes))
    DimPlot(objects)
    
    p1 = DimPlot(objects, cells = bm_cells , pt.size = 0.3, raster = FALSE)
    p2 = DimPlot(objects, cells = cpj_cells, pt.size = 0.3, raster = FALSE)
    p3 = DimPlot(objects, cells = ucb_cells, pt.size = 0.3, raster = FALSE)
    p4 = DimPlot(objects, cells = wj_cells,  pt.size = 0.3, raster = FALSE)
    
    p1 = ggpar(p1, title = "BM-MSC") 
    p2 = ggpar(p2, title = "CPJ-MSC")
    p3 = ggpar(p3, title = "UCB-MSC") 
    p4 = ggpar(p4, title = "WJ-MSC")
    
    ggsave(ggarrange(p1, p2, p3, p4, ncol = 4, labels = "auto"), width = 24, height = 6, dpi = 300,
           device = "tiff", compression = "lzw",
           filename = paste0("./outputs/",project_name,"_subset_cellcycle.tiff"))
} # Use before integration separately
run_seurat_cellcyle_chisq = function(objects, project_name){
    
    tmp_df = as.data.frame(objects[[c("orig.ident","Phase")]])
    tmp_df$Phase = as.character(tmp_df$Phase)
    
    p1 = ggbarstats(tmp_df, x = "orig.ident", y = "Phase", xlab = "Phase", ylab = "Percentage") +
        theme(axis.text.x = element_text(size = 16, color = "black"), axis.text.y = element_text(size = 14, color = "black")) +
        theme(axis.title.x = element_text(size = 18, color = "black"), axis.title.y = element_text(size = 16, color = "black")) +
        scale_fill_brewer(palette="Pastel1")
    
    ggsave(p1, filename = paste0("./outputs/",project_name,"_cellchisq_01.tiff"), width = 12, height = 8, dpi = 600,
           device = "tiff", compression = "lzw")
} # Use before integration separately

#### 4. Lineage ####
run_seurat_lineage_score = function(objects){
    
    Idents(objects) = objects@meta.data$orig.ident
    DefaultAssay(objects) = "SCT"
    
    chondrogenic_markers = list(c("HMGB1", "HMGB2", "DCN", "F3", "MDK"))
    adipogenic_markers = list(c("FTH1", "TAGLN", "FKBP1A", "ACTG2", "TXN"))
    osteogenic_markers = list(c("BGN", "HAPLN1", "FHL1", "VCAN", "GDF15"))
    
    chondrogenic_score = AddModuleScore(objects, chondrogenic_markers)$Cluster1
    adipogenic_score = AddModuleScore(objects, adipogenic_markers)$Cluster1
    osteogenic_score = AddModuleScore(objects, osteogenic_markers)$Cluster1
    
    objects$chondrogenic_strength = chondrogenic_score
    objects$adipogenic_strength = adipogenic_score
    objects$osteogenic_strength = osteogenic_score
    
    return(objects)
}
run_seurat_lineage_assign = function(objects){
    
    chondrogenic_score_95th = quantile(objects$chondrogenic_strength, 0.95)
    adipogenic_score_95th = quantile(objects$adipogenic_strength, 0.95)
    osteogenic_score_95th = quantile(objects$osteogenic_strength, 0.95)
    
    objects$chondrogenic_95th = ifelse(objects$chondrogenic_strength > chondrogenic_score_95th, "CC","NC")
    objects$adipogenic_95th = ifelse(objects$adipogenic_strength > adipogenic_score_95th, "CA","NA")
    objects$osteogenic_95th = ifelse(objects$osteogenic_strength > osteogenic_score_95th, "CO","NO")
    
    return(objects)
}
run_seurat_lineage_chisq = function(objects, project_name){
    
    tmp_df = as.data.frame(objects[[c("orig.ident","chondrogenic_95th")]])
    tmp_df$chondrogenic_95th = as.character(tmp_df$chondrogenic_95th)
    
    p1 = ggbarstats(tmp_df, x = "orig.ident", y = "chondrogenic_95th", xlab = "Chondrogenic lineages", ylab = "Percentage",
                    results.subtitle = FALSE) +
        theme(axis.text.x = element_text(size = 16, color = "black"), axis.text.y = element_text(size = 14, color = "black")) +
        theme(axis.title.x = element_text(size = 18, color = "black"), axis.title.y = element_text(size = 16, color = "black")) +
        scale_fill_brewer(palette="Pastel1")
    
    tmp_df = as.data.frame(objects[[c("orig.ident","adipogenic_95th")]])
    tmp_df$adipogenic_95th = as.character(tmp_df$adipogenic_95th)
    
    p2 = ggbarstats(tmp_df, x = "orig.ident", y = "adipogenic_95th", xlab = "Adipogenic lineages", ylab = "Percentage",
                    results.subtitle = FALSE) +
        theme(axis.text.x = element_text(size = 16, color = "black"), axis.text.y = element_text(size = 14, color = "black")) +
        theme(axis.title.x = element_text(size = 18, color = "black"), axis.title.y = element_text(size = 16, color = "black")) +
        scale_fill_brewer(palette="Pastel1")
    
    tmp_df = as.data.frame(objects[[c("orig.ident","osteogenic_95th")]])
    tmp_df$osteogenic_95th = as.character(tmp_df$osteogenic_95th)
    
    p3 = ggbarstats(tmp_df, x = "orig.ident", y = "osteogenic_95th", xlab = "Osteogenic lineages", ylab = "Percentage",
                    results.subtitle = FALSE) +
        theme(axis.text.x = element_text(size = 16, color = "black"), axis.text.y = element_text(size = 14, color = "black")) +
        theme(axis.title.x = element_text(size = 18, color = "black"), axis.title.y = element_text(size = 16, color = "black")) +
        scale_fill_brewer(palette="Pastel1")
    
    ggsave(p1, filename = paste0("./outputs/",project_name,"_chondrogenic_OR.tiff"), width = 12, height = 8, dpi = 600,
           device = "tiff", compression = "lzw")
    ggsave(p2, filename = paste0("./outputs/",project_name,"_adipogenic_OR.tiff"), width = 12, height = 8, dpi = 600,
           device = "tiff", compression = "lzw")
    ggsave(p3, filename = paste0("./outputs/",project_name,"_osteogenic_OR.tiff"), width = 12, height = 8, dpi = 600,
           device = "tiff", compression = "lzw")
}

#### 5. Clustering ####

run_seurat_umap = function(objects, project_name, res){
    
    objects = RunPCA(objects, verbose = FALSE)
    objects = RunUMAP(objects, dims = 1:30)
    objects = FindNeighbors(objects, dims = 1:30)
    objects = FindClusters(objects, resolution = res)
    
    p1 = DimPlot(objects, reduction = "umap", group.by = "orig.ident", pt.size = 0.3, raster = FALSE)
    p2 = DimPlot(objects, reduction = "umap", group.by = "ident", label = TRUE, pt.size = 0.3, raster = FALSE, label.size = 5)
    p1 = ggpar(p1, title = "") 
    p2 = ggpar(p2, title = "")
    
    ggsave(p1, filename = paste0("./outputs/",project_name,"_umap_origin.tiff"), width = 8, height = 8, dpi = 600,
           device = "tiff", compression = "lzw")
    ggsave(p2, filename = paste0("./outputs/",project_name,"_umap_cluster.tiff"),  width = 8, height = 8, dpi = 600,
           device = "tiff", compression = "lzw")
    ggsave(ggarrange(p1, p2, ncol = 2, labels = "auto"), width = 16, height = 8, dpi = 600,
           device = "tiff", compression = "lzw",
           filename = paste0("./outputs/",project_name,"_umap_origincluster.tiff"))
    
    return(objects)
}
run_seurat_tsne = function(objects, project_name, res){
    
    objects = RunPCA(objects, verbose = FALSE)
    objects = RunTSNE(objects, dims = 1:30)
    objects = FindNeighbors(objects, dims = 1:30)
    objects = FindClusters(objects, resolution = res)
    
    p1 = DimPlot(objects, reduction = "tsne", group.by = "orig.ident", pt.size = 0.3, raster = FALSE)
    p2 = DimPlot(objects, reduction = "tsne", group.by = "ident", label = TRUE, pt.size = 0.3, raster = FALSE, label.size = 5)
    p1 = ggpar(p1, title = "") 
    p2 = ggpar(p2, title = "")
    
    ggsave(p1, filename = paste0("./outputs/",project_name,"_tsne_origin.tiff"), width = 8, height = 8, dpi = 600,
           device = "tiff", compression = "lzw")
    ggsave(p2, filename = paste0("./outputs/",project_name,"_tsne_cluster.tiff"),  width = 8, height = 8, dpi = 600,
           device = "tiff", compression = "lzw")
    ggsave(ggarrange(p1, p2, ncol = 2, labels = "auto"), width = 16, height = 8, dpi = 600,
           device = "tiff", compression = "lzw",
           filename = paste0("./outputs/",project_name,"_tsne_origincluster.tiff"))
    
    return(objects)
}

run_seurat_umap_subset = function(objects, project_name, res){
    
    Idents(objects) = objects@meta.data$orig.ident
    
    bm_cells  = WhichCells(objects, idents = "BM-MSC" )
    cpj_cells = WhichCells(objects, idents = "CPJ-MSC")
    ucb_cells = WhichCells(objects, idents = "UCB-MSC")
    wj_cells  = WhichCells(objects, idents = "WJ-MSC" )
    
    objects = RunPCA(objects, verbose = FALSE)
    objects = RunUMAP(objects, dims = 1:30)
    objects = FindNeighbors(objects, dims = 1:30)
    objects = FindClusters(objects, resolution = res)
    
    p1 = DimPlot(objects, reduction = "umap", group.by = "ident", label = TRUE, cells = bm_cells , pt.size = 0.3, raster = FALSE, label.size = 5) + NoLegend()
    p2 = DimPlot(objects, reduction = "umap", group.by = "ident", label = TRUE, cells = cpj_cells, pt.size = 0.3, raster = FALSE, label.size = 5) + NoLegend()
    p3 = DimPlot(objects, reduction = "umap", group.by = "ident", label = TRUE, cells = ucb_cells, pt.size = 0.3, raster = FALSE, label.size = 5) + NoLegend()
    p4 = DimPlot(objects, reduction = "umap", group.by = "ident", label = TRUE, cells = wj_cells,  pt.size = 0.3, raster = FALSE, label.size = 5) + NoLegend()
    
    p1 = ggpar(p1, title = "BM-MSC") 
    p2 = ggpar(p2, title = "CPJ-MSC")
    p3 = ggpar(p3, title = "UCB-MSC") 
    p4 = ggpar(p4, title = "WJ-MSC")
    
    ggsave(ggarrange(p1, p2, p3, p4, ncol = 4, labels = "auto"), width = 24, height = 6, dpi = 600,
           device = "tiff", compression = "lzw",
           filename = paste0("./outputs/",project_name,"_umap_subset_origincluster.tiff"))
    
    return(objects)
}
run_seurat_tsne_subset = function(objects, project_name, res){
    
    Idents(objects) = objects@meta.data$orig.ident
    
    bm_cells  = WhichCells(objects, idents = "BM-MSC" )
    cpj_cells = WhichCells(objects, idents = "CPJ-MSC")
    ucb_cells = WhichCells(objects, idents = "UCB-MSC")
    wj_cells  = WhichCells(objects, idents = "WJ-MSC" )
    
    objects = RunPCA(objects, verbose = FALSE)
    objects = RunTSNE(objects, dims = 1:30)
    objects = FindNeighbors(objects, dims = 1:30)
    objects = FindClusters(objects, resolution = res)
    
    p1 = DimPlot(objects, reduction = "tsne", group.by = "ident", label = TRUE, cells = bm_cells , pt.size = 0.3, raster = FALSE, label.size = 5) + NoLegend()
    p2 = DimPlot(objects, reduction = "tsne", group.by = "ident", label = TRUE, cells = cpj_cells, pt.size = 0.3, raster = FALSE, label.size = 5) + NoLegend()
    p3 = DimPlot(objects, reduction = "tsne", group.by = "ident", label = TRUE, cells = ucb_cells, pt.size = 0.3, raster = FALSE, label.size = 5) + NoLegend()
    p4 = DimPlot(objects, reduction = "tsne", group.by = "ident", label = TRUE, cells = wj_cells,  pt.size = 0.3, raster = FALSE, label.size = 5) + NoLegend()
    
    p1 = ggpar(p1, title = "BM-MSC") 
    p2 = ggpar(p2, title = "CPJ-MSC")
    p3 = ggpar(p3, title = "UCB-MSC") 
    p4 = ggpar(p4, title = "WJ-MSC")
    
    ggsave(ggarrange(p1, p2, p3, p4, ncol = 4, labels = "auto"), width = 24, height = 6, dpi = 600,
           device = "tiff", compression = "lzw",
           filename = paste0("./outputs/",project_name,"_tsne_subset_origincluster.tiff"))
    
    return(objects)
}

run_seurat_umap_lineage = function(objects, project_name, res){
    
    DefaultAssay(objects) = "integrated"
    objects = RunPCA(objects, verbose = FALSE)
    objects = RunUMAP(objects, dims = 1:30)
    objects = FindNeighbors(objects, dims = 1:30)
    objects = FindClusters(objects, resolution = res)
    
    p1 = FeaturePlot(objects, features = c("chondrogenic_strength", "adipogenic_strength", "osteogenic_strength"),
                     col = c("blue","grey","red"), ncol = 3, pt.size = 0.4)
    
    ggsave(p1, filename = paste0("./outputs/",project_name,"_umap_lineage.tiff"), width = 24, height = 8, dpi = 300,
           device = "tiff", compression = "lzw")

}
run_seurat_tsne_lineage = function(objects, project_name, res){
    
    DefaultAssay(objects) = "integrated"
    objects = RunPCA(objects, verbose = FALSE)
    objects = RunTSNE(objects, dims = 1:30)
    objects = FindNeighbors(objects, dims = 1:30)
    objects = FindClusters(objects, resolution = res)
    
    p1 = FeaturePlot(objects, reduction = "tsne", features = c("chondrogenic_strength", "adipogenic_strength", "osteogenic_strength"),
                     col = c("blue","grey","red"), ncol = 3, pt.size = 0.4)
    
    ggsave(p1, filename = paste0("./outputs/",project_name,"_tsne_lineage.tiff"), width = 24, height = 8, dpi = 300,
           device = "tiff", compression = "lzw")
    
}

run_seurat_cluster_source_chisq = function(objects, project_name){
    
    tmp_df = as.data.frame(objects[[c("orig.ident","integrated_snn_res.0.1", "integrated_snn_res.0.6")]])
    tmp_df$integrated_snn_res.0.1 = as.integer(as.character(tmp_df$integrated_snn_res.0.1))
    tmp_df$integrated_snn_res.0.6 = as.integer(as.character(tmp_df$integrated_snn_res.0.6))
    
    p1 = ggbarstats(tmp_df, x = "orig.ident", y = "integrated_snn_res.0.1", xlab = "Clusters", ylab = "Percentage",
                    results.subtitle = FALSE) +
        theme(axis.text.x = element_text(size = 16, color = "black"), axis.text.y = element_text(size = 14, color = "black")) +
        theme(axis.title.x = element_text(size = 18, color = "black"), axis.title.y = element_text(size = 16, color = "black")) +
        scale_fill_brewer(palette="Pastel1")
    
    p2 = ggbarstats(tmp_df, x = "orig.ident", y = "integrated_snn_res.0.6", xlab = "Clusters", ylab = "Percentage",
                    results.subtitle = FALSE) +
        theme(axis.text.x = element_text(size = 16, color = "black"), axis.text.y = element_text(size = 14, color = "black")) +
        theme(axis.title.x = element_text(size = 18, color = "black"), axis.title.y = element_text(size = 16, color = "black")) +
        scale_fill_brewer(palette="Pastel1")
    
    ggsave(p1, filename = paste0("./outputs/",project_name,"_chisq_01.tiff"), width = 12, height = 8, dpi = 600,
           device = "tiff", compression = "lzw")
    ggsave(p2, filename = paste0("./outputs/",project_name,"_chisq_06.tiff"), width = 16, height = 8, dpi = 600,
           device = "tiff", compression = "lzw")
}
run_seurat_cluster_lineage_chondrogenic_chisq = function(objects, project_name){
    
    tmp_df = as.data.frame(objects[[c("chondrogenic_95th","integrated_snn_res.0.1", "integrated_snn_res.0.6")]])
    tmp_df$integrated_snn_res.0.1 = as.integer(as.character(tmp_df$integrated_snn_res.0.1))
    tmp_df$integrated_snn_res.0.6 = as.integer(as.character(tmp_df$integrated_snn_res.0.6))
    
    p1 = ggbarstats(tmp_df, x = "chondrogenic_95th", y = "integrated_snn_res.0.1", xlab = "Chondrongenic lineage", ylab = "Percentage"
                    ,results.subtitle = FALSE) +
        theme(axis.text.x = element_text(size = 16, color = "black"), axis.text.y = element_text(size = 14, color = "black")) +
        theme(axis.title.x = element_text(size = 18, color = "black"), axis.title.y = element_text(size = 16, color = "black")) +
        scale_fill_brewer(palette="Pastel1")
    
    p2 = ggbarstats(tmp_df, x = "chondrogenic_95th", y = "integrated_snn_res.0.6", xlab = "Chondrongenic lineage", ylab = "Percentage",
                    results.subtitle = FALSE) +
        theme(axis.text.x = element_text(size = 16, color = "black"), axis.text.y = element_text(size = 14, color = "black")) +
        theme(axis.title.x = element_text(size = 18, color = "black"), axis.title.y = element_text(size = 16, color = "black")) +
        scale_fill_brewer(palette="Pastel1")
    
    ggsave(p1, filename = paste0("./outputs/",project_name,"_chondrogenic_chisq_01.tiff"), width = 12, height = 8, dpi = 600,
           device = "tiff", compression = "lzw")
    ggsave(p2, filename = paste0("./outputs/",project_name,"_chondrogenic_chisq_06.tiff"), width = 16, height = 8, dpi = 600,
           device = "tiff", compression = "lzw")
}
run_seurat_cluster_lineage_adipogenic_chisq = function(objects, project_name){
    
    tmp_df = as.data.frame(objects[[c("adipogenic_95th","integrated_snn_res.0.1", "integrated_snn_res.0.6")]])
    tmp_df$integrated_snn_res.0.1 = as.integer(as.character(tmp_df$integrated_snn_res.0.1))
    tmp_df$integrated_snn_res.0.6 = as.integer(as.character(tmp_df$integrated_snn_res.0.6))
    
    p1 = ggbarstats(tmp_df, x = "adipogenic_95th", y = "integrated_snn_res.0.1", xlab = "Adipogenic lineage", ylab = "Percentage",
                    results.subtitle = FALSE) +
        theme(axis.text.x = element_text(size = 16, color = "black"), axis.text.y = element_text(size = 14, color = "black")) +
        theme(axis.title.x = element_text(size = 18, color = "black"), axis.title.y = element_text(size = 16, color = "black")) +
        scale_fill_brewer(palette="Pastel1")
    
    p2 = ggbarstats(tmp_df, x = "adipogenic_95th", y = "integrated_snn_res.0.6", xlab = "Adipogenic lineage", ylab = "Percentage",
                    results.subtitle = FALSE) +
        theme(axis.text.x = element_text(size = 16, color = "black"), axis.text.y = element_text(size = 14, color = "black")) +
        theme(axis.title.x = element_text(size = 18, color = "black"), axis.title.y = element_text(size = 16, color = "black")) +
        scale_fill_brewer(palette="Pastel1")
    
    ggsave(p1, filename = paste0("./outputs/",project_name,"_adipogenic_chisq_01.tiff"), width = 12, height = 8, dpi = 600,
           device = "tiff", compression = "lzw")
    ggsave(p2, filename = paste0("./outputs/",project_name,"_adipogenic_chisq_06.tiff"), width = 16, height = 8, dpi = 600,
           device = "tiff", compression = "lzw")
}
run_seurat_cluster_lineage_osteogenic_chisq = function(objects, project_name){
    
    tmp_df = as.data.frame(objects[[c("osteogenic_95th","integrated_snn_res.0.1", "integrated_snn_res.0.6")]])
    tmp_df$integrated_snn_res.0.1 = as.integer(as.character(tmp_df$integrated_snn_res.0.1))
    tmp_df$integrated_snn_res.0.6 = as.integer(as.character(tmp_df$integrated_snn_res.0.6))
    
    p1 = ggbarstats(tmp_df, x = "osteogenic_95th", y = "integrated_snn_res.0.1", xlab = "Osteogenic lineage", ylab = "Percentage",
                    results.subtitle = FALSE) +
        theme(axis.text.x = element_text(size = 16, color = "black"), axis.text.y = element_text(size = 14, color = "black")) +
        theme(axis.title.x = element_text(size = 18, color = "black"), axis.title.y = element_text(size = 16, color = "black")) +
        scale_fill_brewer(palette="Pastel1")
    
    p2 = ggbarstats(tmp_df, x = "osteogenic_95th", y = "integrated_snn_res.0.6", xlab = "Osteogenic lineage", ylab = "Percentage",
                    results.subtitle = FALSE) +
        theme(axis.text.x = element_text(size = 16, color = "black"), axis.text.y = element_text(size = 14, color = "black")) +
        theme(axis.title.x = element_text(size = 18, color = "black"), axis.title.y = element_text(size = 16, color = "black")) +
        scale_fill_brewer(palette="Pastel1")
    
    ggsave(p1, filename = paste0("./outputs/",project_name,"_osteogenic_chisq_01.tiff"), width = 12, height = 8, dpi = 600,
           device = "tiff", compression = "lzw")
    ggsave(p2, filename = paste0("./outputs/",project_name,"_osteogenic_chisq_06.tiff"), width = 16, height = 8, dpi = 600,
           device = "tiff", compression = "lzw")
}

#### 6. MSC effector scores ####

run_seurat_ALSMSC_score = function(objects){
    
    Idents(objects) = objects@meta.data$orig.ident
    DefaultAssay(objects) = "SCT"

    als_pos_markers = list(c("TGFB1","TGFB2","TGFB3","IL6"))
    als_neg_markers = list(c("CCL2","TNF"))
    msc_markers = list(c("TGFB1", "VEGFA", "IDO1", "IL6", "CXCL8", "PTGES2", "CXCL12", "HGF", "CD274", "TNFAIP6", "IL1RN"))
    #msc_markers = list(c("TGFB1", "VEGFA", "CXCL8", "PTGES2", "CXCL12", "HGF", "TNFAIP6"))
    
    als_pos_score = AddModuleScore(objects, als_pos_markers)$Cluster1
    als_neg_score = AddModuleScore(objects, als_neg_markers)$Cluster1
    msc_score = AddModuleScore(objects, msc_markers)$Cluster1
    als_score = als_pos_score - als_neg_score
    
    objects$als_effector_strength = als_score
    objects$msc_effector_strength = msc_score
    
    return(objects)
}
run_seurat_umap_effector_score = function(objects, project_name, res){
    
    DefaultAssay(objects) = "integrated"
    objects = RunPCA(objects, verbose = FALSE)
    objects = RunUMAP(objects, dims = 1:30)
    objects = FindNeighbors(objects, dims = 1:30)
    objects = FindClusters(objects, resolution = res)
    
    p1 = FeaturePlot(objects, features = c("als_effector_strength", "msc_effector_strength"),
                     col = c("blue","grey","red"), ncol = 2, pt.size = 0.4)
    
    ggsave(p1, filename = paste0("./outputs/",project_name,"_umap_effector_strength.tiff"), width = 16, height = 8, dpi = 600,
           device = "tiff", compression = "lzw")
    
}
run_seurat_tsne_effector_score = function(objects, project_name, res){
    
    DefaultAssay(objects) = "integrated"
    objects = RunPCA(objects, verbose = FALSE)
    objects = RunTSNE(objects, dims = 1:30)
    objects = FindNeighbors(objects, dims = 1:30)
    objects = FindClusters(objects, resolution = res)
    
    p1 = FeaturePlot(objects, reduction = "tsne", features = c("als_effector_strength", "msc_effector_strength"),
                     col = c("blue","grey","red"), ncol = 2, pt.size = 0.4)
    
    ggsave(p1, filename = paste0("./outputs/",project_name,"_tsne_effector_strength.tiff"), width = 16, height = 8, dpi = 600,
           device = "tiff", compression = "lzw")
    
}
run_seurat_msc_score_boxplot = function(objects, project_name){
    
    tmp_df = as.data.frame(objects[[c("msc_effector_strength","osteogenic_strength", "adipogenic_strength", "chondrogenic_strength", "orig.ident")]])
    
    p1 = ggbetweenstats(tmp_df, y = "msc_effector_strength", x = "orig.ident", xlab = "Origin", ylab = "MSC effector strength",
                        color = "orig.ident", pairwise.comparisons = FALSE, results.subtitle = FALSE, type = "p", tr = 0, k = 4) +
        theme(axis.text.x = element_text(size = 16, color = "black"), axis.text.y = element_text(size = 14, color = "black")) +
        theme(axis.title.x = element_text(size = 18, color = "black"), axis.title.y = element_text(size = 16, color = "black")) +
        scale_fill_brewer(palette="Pastel1")
    
    ggsave(p1, filename = paste0("./outputs/",project_name,"_origin_msc_strength.tiff"), width = 16, height = 8, dpi = 300,
           device = "tiff", compression = "lzw")
    
    tmp_df = as.data.frame(objects[[c("msc_effector_strength","integrated_snn_res.0.6")]])
    tmp_df$integrated_snn_res.0.6 = as.integer(as.character(tmp_df$integrated_snn_res.0.6))
    
    p2 = ggbetweenstats(tmp_df, y = "msc_effector_strength", x = "integrated_snn_res.0.6", xlab = "Clusters", ylab = "MSC effector strength",
                       color = "integrated_snn_res.0.6", pairwise.comparisons = FALSE, results.subtitle = FALSE, type = "p", tr = 0, k = 4) +
        theme(axis.text.x = element_text(size = 16, color = "black"), axis.text.y = element_text(size = 14, color = "black")) +
        theme(axis.title.x = element_text(size = 18, color = "black"), axis.title.y = element_text(size = 16, color = "black")) +
        scale_color_manual(values = colorRampPalette(brewer.pal(8, "Set2"))(20)) + 
        scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Set2"))(20))
    
    ggsave(p2, filename = paste0("./outputs/",project_name,"_HC_msc_strength.tiff"), width = 24, height = 12, dpi = 600,
           device = "tiff", compression = "lzw")
    
    tmp_df = as.data.frame(objects[[c("msc_effector_strength","integrated_snn_res.0.1")]])
    tmp_df$integrated_snn_res.0.1 = as.integer(as.character(tmp_df$integrated_snn_res.0.1))
    
    p3 = ggbetweenstats(tmp_df, y = "msc_effector_strength", x = "integrated_snn_res.0.1", xlab = "Clusters", ylab = "MSC effector strength",
                       color = "integrated_snn_res.0.1", pairwise.comparisons = FALSE, results.subtitle = FALSE, type = "p", tr = 0, k = 4) +
        theme(axis.text.x = element_text(size = 16, color = "black"), axis.text.y = element_text(size = 14, color = "black")) +
        theme(axis.title.x = element_text(size = 18, color = "black"), axis.title.y = element_text(size = 16, color = "black")) +
        scale_color_manual(values = colorRampPalette(brewer.pal(8, "Set2"))(20)) + 
        scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Set2"))(20))
    
    ggsave(p3, filename = paste0("./outputs/",project_name,"_LC_msc_strength.tiff"), width = 18, height = 12, dpi = 600,
           device = "tiff", compression = "lzw")
    
    tmp_df = as.data.frame(objects[[c("msc_effector_strength", "osteogenic_strength", "adipogenic_strength", "chondrogenic_strength", "orig.ident")]])
    
    p4 = ggbetweenstats(tmp_df, y = "chondrogenic_strength", x = "orig.ident", xlab = "Origin", ylab = "Chondrogenic strength",
                        color = "orig.ident", pairwise.comparisons = FALSE, results.subtitle = FALSE, type = "p", tr = 0, k = 4) +
        theme(axis.text.x = element_text(size = 16, color = "black"), axis.text.y = element_text(size = 14, color = "black")) +
        theme(axis.title.x = element_text(size = 18, color = "black"), axis.title.y = element_text(size = 16, color = "black")) +
        scale_fill_brewer(palette="Pastel1")
    
    ggsave(p4, filename = paste0("./outputs/",project_name,"_origin_chondrogenic_strength.tiff"), width = 16, height = 8, dpi = 300,
           device = "tiff", compression = "lzw")
    
    p5 = ggbetweenstats(tmp_df, y = "osteogenic_strength", x = "orig.ident", xlab = "Origin", ylab = "Osteogenic strength",
                        color = "orig.ident", pairwise.comparisons = FALSE, results.subtitle = FALSE, type = "p", tr = 0, k = 4) +
        theme(axis.text.x = element_text(size = 16, color = "black"), axis.text.y = element_text(size = 14, color = "black")) +
        theme(axis.title.x = element_text(size = 18, color = "black"), axis.title.y = element_text(size = 16, color = "black")) +
        scale_fill_brewer(palette="Pastel1")
    
    ggsave(p5, filename = paste0("./outputs/",project_name,"_origin_osteogenic_strength.tiff"), width = 16, height = 8, dpi = 300,
           device = "tiff", compression = "lzw")
    
    p6 = ggbetweenstats(tmp_df, y = "adipogenic_strength", x = "orig.ident", xlab = "Origin", ylab = "Adipogenic strength",
                        color = "orig.ident", pairwise.comparisons = FALSE, results.subtitle = FALSE, type = "p", tr = 0, k = 4) +
        theme(axis.text.x = element_text(size = 16, color = "black"), axis.text.y = element_text(size = 14, color = "black")) +
        theme(axis.title.x = element_text(size = 18, color = "black"), axis.title.y = element_text(size = 16, color = "black")) +
        scale_fill_brewer(palette="Pastel1")
    
    ggsave(p6, filename = paste0("./outputs/",project_name,"_origin_adipogenic_strength.tiff"), width = 16, height = 8, dpi = 300,
           device = "tiff", compression = "lzw")
}

#### 7. Markers ####
run_seurat_normalization = function(objects){
    DefaultAssay(objects) = "RNA"
    objects = NormalizeData(objects)
    objects = FindVariableFeatures(objects, selection.method = "vst", nfeatures = 2000)
    all_genes = rownames(objects)
    objects = ScaleData(objects, features = all_genes)
    return(objects)
}
run_seurat_ridge_markers = function(objects, project_name){
    Idents(objects) = objects@meta.data$orig.ident
    msc_markers = c("ITGB1", "CD44", "NT5E", "THY1", "ENG", "CD34", "PTPRC", "HLA-DRB1")
    msc_protein = c("CD29", "CD44", "CD73", "CD90", "CD105", "CD34", "CD45", "HLA-DRB1")
    
    for (i in 1:length(msc_markers)){
        assign(x = paste0("tmp",i),value = (RidgePlot(objects, features = msc_markers[i], cols = palette("pastel1")) +
                                                theme(axis.title.y = element_blank()) +
                                                theme(legend.position = "none") +
                                                labs(title = msc_protein[i])), envir = .GlobalEnv)
    }
    all_plot = ggarrange(tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, ncol = 5, nrow = 2, labels = "auto")
    ggsave(all_plot, width = 16, height = 8, dpi = 300, device = "tiff",
           filename = paste0("./outputs/",project_name,"_ridge_markers.tiff"))
}
run_seurat_violin_markers = function(objects, project_name){
    Idents(objects) = objects@meta.data$orig.ident
    msc_markers = c("ITGB1", "CD44", "NT5E", "THY1", "ENG", "CD34", "PTPRC", "HLA-DRB1")
    msc_protein = c("CD29", "CD44", "CD73", "CD90", "CD105", "CD34", "CD45", "HLA-DRB1")
    
    for (i in 1:length(msc_markers)){
        assign(x = paste0("tmp",i),value = (VlnPlot(objects, features = msc_markers[i], cols = palette("pastel1")) +
                                                theme(axis.title.y = element_blank()) +
                                                theme(legend.position = "none") +
                                                labs(title = msc_protein[i])), envir = .GlobalEnv)
    }
    all_plot = ggarrange(tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, ncol = 5, nrow = 2, labels = "auto")
    ggsave(all_plot, width = 16, height = 8, dpi = 300, device = "tiff",
           filename = paste0("./outputs/",project_name,"_violin_markers.tiff"))
}

run_seurat_violin_ALS_markers = function(objects, project_name){
    Idents(objects) = objects@meta.data$orig.ident
    msc_markers = c("TGFB1","TGFB2","TGFB3","IL6","CCL2","TNF")
    msc_protein = c("TGF-B1","TGF-B2","TGF-B3","IL-6","MCP-1","TNF-a")
    
    for (i in 1:length(msc_markers)){
        assign(x = paste0("tmp",i),value = (VlnPlot(objects, features = msc_markers[i], cols = palette("pastel1")) +
                                                theme(axis.title.y = element_blank()) +
                                                theme(legend.position = "none") +
                                                stat_summary(fun = mean, geom = "point", color = "red", size = 3) +
                                                labs(title = msc_protein[i])), envir = .GlobalEnv)
    }
    all_plot = ggarrange(tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, ncol = 3, nrow = 2, labels = "auto")
    ggsave(all_plot, width = 16, height = 10, dpi = 300, device = "tiff",
           filename = paste0("./outputs/",project_name,"_violin_ALS_markers.tiff"))
}
run_seurat_violin_MSC_markers = function(objects, project_name){
    Idents(objects) = objects@meta.data$orig.ident
    msc_markers = c("TGFB1", "VEGFA", "IDO1", "IL6", "CXCL8", "PTGES2", "CXCL12", "HGF", "CD274", "TNFAIP6", "IL1RN")
    msc_protein = c("TGF-B1", "VEGF-A", "IDO-1", "IL-6", "IL-8", "PGE2", "SDF-1", "HGF", "PD-L1", "TSF-6", "IL-1Ra")
    
    for (i in 1:length(msc_markers)){
        assign(x = paste0("tmp",i),value = (VlnPlot(objects, features = msc_markers[i], cols = palette("pastel1")) +
                                                theme(axis.title.y = element_blank()) +
                                                theme(legend.position = "none") +
                                                stat_summary(fun = mean, geom = "point", color = "red", size = 3) +
                                                labs(title = msc_protein[i])), envir = .GlobalEnv)
    }
    all_plot = ggarrange(tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11,
                         ncol = 3, nrow = 4, labels = "auto")
    ggsave(all_plot, width = 16, height = 20, dpi = 300, device = "tiff",
           filename = paste0("./outputs/",project_name,"_violin_MSC_markers.tiff"))
}
run_seurat_violin_ALSMSC_score = function(objects, project_name){
    
    Idents(objects) = objects@meta.data$orig.ident
    
    als_pos_markers = list(c("TGFB1","TGFB2","TGFB3","IL6"))
    als_neg_markers = list(c("CCL2","TNF"))
    msc_markers = list(c("TGFB1", "VEGFA", "IDO1", "IL6", "CXCL8", "PTGES2", "CXCL12", "HGF", "CD274", "TNFAIP6", "IL1RN"))
    
    als_pos_score = AddModuleScore(objects, als_pos_markers)$Cluster1
    als_neg_score = AddModuleScore(objects, als_neg_markers)$Cluster1
    msc_score = AddModuleScore(objects, msc_markers)$Cluster1
    als_score = als_pos_score - als_neg_score
    
    objects$Cluster1 = msc_score
    
    p1 = VlnPlot(objects, features = "Cluster1", cols = palette("pastel1")) +
        theme(axis.title.y = element_blank()) +
        theme(legend.position = "none") +
        stat_summary(fun = mean, geom = "point", color = "red", size = 3) +
        labs(title = "MSC effect score")
    
    objects$Cluster1 = als_score
    
    
    p2 = VlnPlot(objects, features = "Cluster1", cols = palette("pastel1")) +
        theme(axis.title.y = element_blank()) +
        theme(legend.position = "none") +
        stat_summary(fun = mean, geom = "point", color = "red", size = 3) +
        labs(title = "ALS effect score")
    
    ggsave(p1, width = 10, height = 8, dpi = 600, device = "tiff", compression = "lzw",
           filename = paste0("./outputs/",project_name,"_violin_MSC_effect.tiff"))
    ggsave(p2, width = 10, height = 8, dpi = 600, device = "tiff", compression = "lzw",
           filename = paste0("./outputs/",project_name,"_violin_ALS_effect.tiff"))
}

run_seurat_markers = function(objects, project_name){
    
    Idents(objects) = objects@meta.data$integrated_snn_res.0.1
    tmpres01 = FindAllMarkers(object = objects, assay = "SCT")
    
    Idents(objects) = objects@meta.data$integrated_snn_res.0.6
    tmpres06 = FindAllMarkers(object = objects, assay = "SCT")
    
    Idents(objects) = objects@meta.data$orig.ident
    tmporigin = FindAllMarkers(object = objects, assay = "SCT")
    
    write.xlsx(x = tmpres01, file = paste0(project_name,".xlsx"),
               sheetName = paste0(project_name,"res01"), append = TRUE, showNA = F)
    write.xlsx(x = tmpres06, file = paste0(project_name,".xlsx"),
               sheetName = paste0(project_name,"res06"), append = TRUE, showNA = F)
    write.xlsx(x = tmporigin, file = paste0(project_name,".xlsx"),
               sheetName = paste0(project_name,"origin"), append = TRUE, showNA = F)
    
    tmp_list = list(tmpres01,tmpres06,tmporigin)
    return(tmp_list)
}
run_seurat_markers_all = function(objects, project_name){
    
    #Idents(objects) = objects@meta.data$integrated_snn_res.0.1
    #tmpres01 = FindAllMarkers(object = objects, assay = "RNA", logfc.threshold = 0, min.pct = 0, only.pos = F, test.use = "MAST")
    
    #Idents(objects) = objects@meta.data$integrated_snn_res.0.6
    #tmpres06 = FindAllMarkers(object = objects, assay = "RNA", logfc.threshold = 0, min.pct = 0, only.pos = F, test.use = "MAST")
    
    Idents(objects) = objects@meta.data$orig.ident
    tmporigin = FindAllMarkers(object = objects, assay = "RNA", logfc.threshold = 0, min.pct = 0, only.pos = F, test.use = "MAST")
    
    #write.xlsx(x = tmpres01, file = paste0(project_name,"_allgenes_LC.xlsx"),
    #           sheetName = paste0(project_name,"res01"), append = TRUE, showNA = F)
    #write.xlsx(x = tmpres06, file = paste0(project_name,"_allgenes_HC.xlsx"),
    #           sheetName = paste0(project_name,"res06"), append = TRUE, showNA = F)
    write.xlsx(x = tmporigin, file = paste0(project_name,"_allgenes_ORIGIN.xlsx"),
               sheetName = paste0(project_name,"origin"), append = TRUE, showNA = F)
    
    #tmp_list = list(tmpres01,tmpres06,tmporigin)
    #return(tmpres01)
}





#### 8. DEG functions ####
run_DEG_profiler = function(DEGs, sheets, project_name){
    
    tmp = read.xlsx(file = DEGs, sheetIndex = sheets)
    tmp$NA. = NULL
    tmp$ENTREZID = select(org.Hs.eg.db, keys=tmp$gene, columns="ENTREZID", keytype="SYMBOL")$ENTREZID
    tmp = tmp[complete.cases(tmp),]
    
    for (i in 1:nrow(tmp)){
        if(nchar(tmp[i,6]) == 1){
            tmp[i,6] = paste0("0",tmp[i,6])
        }
    }
    
    clusters = length(unique(tmp$cluster))
        
    for (j in 1:clusters){
        selected_cluster = unique(tmp$cluster)[j]
        sub_tmp = subset(tmp, cluster == selected_cluster)
        sub_tmp_list = sub_tmp[,c("ENTREZID", "avg_log2FC")]
        
        sub_geneList = sub_tmp_list[,2]
        names(sub_geneList) = as.character(sub_tmp_list[,1])
        sub_geneList = sort(sub_geneList, decreasing = TRUE)
        
        sub_ego = enrichGO(gene      = names(sub_geneList),
                       OrgDb         = org.Hs.eg.db,
                       ont           = "MF",
                       pAdjustMethod = "none",
                       qvalueCutoff  = 0.1,
                       readable      = TRUE)
        
        sub_tmp_go = dotplot(sub_ego, showCategory = 20)
        
        ggsave(filename = paste0(project_name,"_GO_", eval(selected_cluster), ".tiff"),
               sub_tmp_go, width = 12, height = 12, device = "tiff", dpi = 300)
        
        sub_kegg = enrichKEGG(gene = names(sub_geneList),
                   organism     = 'hsa',
                   pvalueCutoff = 0.05)
        
        sub_tmp_kegg = dotplot(sub_kegg, showCategory = 20)
        
        ggsave(filename = paste0(project_name,"_KEGG_", eval(selected_cluster), ".tiff"),
               sub_tmp_kegg, width = 12, height = 12, device = "tiff", dpi = 300)
        
    }

}
run_DEG_multiprofiler = function(DEGs, sheets, project_name){
    
    tmp = read.xlsx(file = DEGs, sheetIndex = sheets)
    tmp$NA. = NULL
    colnames(tmp)[5] = "p.adjust"
    tmp$ENTREZID = select(org.Hs.eg.db, keys=tmp$gene, columns="ENTREZID", keytype="SYMBOL")$ENTREZID
    tmp = tmp[complete.cases(tmp),]
    
    for (i in 1:nrow(tmp)){
        if(nchar(tmp[i,6]) == 1){
            tmp[i,6] = paste0("0",tmp[i,6])
        }
    }
    
    p1 = compareCluster(ENTREZID~cluster, data=tmp, fun='enrichGO', OrgDb='org.Hs.eg.db', ont = "MF")
    p2 = compareCluster(ENTREZID~cluster, data=tmp, fun='enrichKEGG')

    p1 = dotplot(p1, showCategory = 5)
    p2 = dotplot(p2, showCategory = 10)
    
    p1
    p2
    
    ggsave(filename = paste0(project_name,"_EGO_MULTI.tiff"),
           p1, width = 12, height = 16, device = "tiff", dpi = 300)
    ggsave(filename = paste0(project_name,"_EKEGG_MULTI.tiff"),
           p2, width = 12, height = 16, device = "tiff", dpi = 300)
}
