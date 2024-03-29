library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
library(ggplot2)
library(patchwork)

counts <- Read10X_h5("~/Data/pbmc_unsorted_3k/pbmc_unsorted_3k_filtered_feature_bc_matrix.h5")
fragments <- "~/Data/pbmc_unsorted_3k/pbmc_unsorted_3k_atac_fragments.tsv.gz"

annotations <- GetGRangesFromEnsDb(EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "hg38"

# call peaks using MACS2
peaks <- CallPeaks(object = fragments)
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")

# remove peaks in blacklist regions
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

# remove low-confidence peaks
peaks <- peaks[peaks$neg_log10qvalue_summit > 5]

# write peaks to file
write.table(x = as.data.frame(peaks), file = "~/Data/pbmc_unsorted_3k/pbmc_multiomic_peaks.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# create object
pbmc <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  project = "coassay",
  assay = "RNA"
)

# quantify peaks
frags <- CreateFragmentObject(
  path = fragments,
  cells = colnames(pbmc)
)
peakcounts <- FeatureMatrix(
  fragments = frags,
  features = peaks,
  cells = colnames(pbmc)
)

pbmc[["ATAC"]] <- CreateChromatinAssay(
  counts = peakcounts,
  min.cells = 5,
  genome = "hg38",
  fragments = frags,
  annotation = annotations
)

pbmc <- subset(
  x = pbmc,
  subset = nCount_ATAC < 100000 &
    nCount_RNA < 25000 &
    nCount_ATAC > 1000 &
    nCount_RNA > 1000
)

DefaultAssay(pbmc) <- "ATAC"
pbmc <- NucleosomeSignal(pbmc)
pbmc <- TSSEnrichment(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 10)
pbmc <- RunTFIDF(pbmc)
pbmc <- RunSVD(pbmc)
pbmc <- RunUMAP(pbmc, reduction = "lsi", dims = 2:40, reduction.name = "umap.atac")

DefaultAssay(pbmc) <- "RNA"
pbmc <- NormalizeData(pbmc)
pbmc <- SCTransform(pbmc, ncells = 5000)
pbmc <- RunPCA(pbmc)
pbmc <- RunUMAP(pbmc, dims = 1:40, reduction.name = "umap.rna")

## Label Transfer using the new pbmc multimodal citeseq dataset as reference dataset

remotes::install_github("jlmelville/uwot")
remotes::install_github("mojaveazure/seurat-disk")
remotes::install_github("satijalab/seurat�, ref = �release/4.0.0") 

library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(patchwork)

reference <- LoadH5Seurat("/home/madads/Data/pbmc_multimodal.h5seurat")

anchors <- FindTransferAnchors(
  reference = reference,
  query = pbmc,
  query.assay= "SCT",
  normalization.method = "SCT",
  reference.reduction = "spca",
  dims = 1:50
)

pbmc <- MapQuery(
  anchorset = anchors,
  query = pbmc,
  reference = reference,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2",
    predicted_ADT = "ADT"
  ),
  reference.reduction = "spca", 
  reduction.model = "wnn.umap"
)

p1 = DimPlot(pbmc, reduction = "ref.umap", group.by = "predicted.celltype.l1", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
p2 = DimPlot(pbmc, reduction = "ref.umap", group.by = "predicted.celltype.l2", label = TRUE, label.size = 3 ,repel = TRUE) + NoLegend()
p1 + p2


saveRDS(object = pbmc, file = "~/Data/pbmc_unsorted_3k/pbmc_unsorted_3k.rds")
