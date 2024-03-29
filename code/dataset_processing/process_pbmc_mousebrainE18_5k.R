library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
library(ggplot2)
library(patchwork)

counts <- Read10X_h5("~/Data/pbmc_mousebrainE18/e18_mouse_brain_fresh_5k_filtered_feature_bc_matrix.h5")
fragments <- "~/Data/pbmc_mousebrainE18/e18_mouse_brain_fresh_5k_atac_fragments.tsv.gz"

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
write.table(x = as.data.frame(peaks), file = "~/Data/pbmc_mousebrainE18/pbmc_mousebrainE18_5k_peaks.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

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

### Label Transfer using allen brain scrna-seq as reference 

# Load in the dataset 

allen_brain<-  readRDS("/home/madads/Data/allen_brain.rds")

transfer.anchors <- FindTransferAnchors(
  reference = allen_brain,
  query = pbmc, reference.assay = "SCT", query.assay = "SCT",
  reduction = 'cca'
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = allen_brain$,
  weight.reduction = pbmc[['pca']],
  dims = 1:30
)

pbmc<- AddMetaData(object = coassay, metadata = predicted.labels)

pbmc[["predicted.id"]]


saveRDS(object = pbmc, file = "~/Data/pbmc_mousebrainE18/pbmc_mousebrainE18_5k.rds")
