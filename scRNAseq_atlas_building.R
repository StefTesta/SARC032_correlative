##### Single-Cell RNA sequencing Workflow for CIBERSORTx signaure building -----
library(Seurat)
library(harmony)
library(dplyr)
library(clustree)
library(sctransform)
library(glmGamPoi)
library(biomaRt)
library(tidyr)
library(limma)
library(edgeR)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(Orthology.eg.db)

set.seed(123)

#### ---- Gather dissociation induced genes (DIGs), TCR genes and Ribosomal protein genes ---- ####
# Lis of DIGs from Van Den Brink et al. -----
mouse_genes <- c("Actg1", "Btg1", "Cxcl1", "Dnajb4", "Errfi1", "H3f3b", "Hspb1", "Irf1", "Klf6", "Mir22hg", 
                 "Nfkbia", "Pcf11", "Pxdc1", "Sdc4", "Srf", "Tpm3", "Usp2", "Gadd45g", "Ankrd1", "Btg2", 
                 "Cyr61", "Dusp1", "Fam132b", "Hipk3", "Hsph1", "Irf8", "Klf9", "Mt1", "Nfkbiz", "Pde4b", 
                 "Rap1b", "Serpine1", "Srsf5", "Tppp3", "Wac", "Hspe1", "Arid5a", "Ccnl1", "Dusp8", "Fos", 
                 "Hsp90aa1", "Id3", "Itpkc", "Litaf", "Mt2", "Nop58", "Per1", "Rassf1", "Skil", "Srsf7", 
                 "Tra2a", "Zc3h12a", "Ier5", "Atf3", "Ccrn4l", "Ddx3x", "Egr1", "Fosb", "Hsp90ab1", "Idi1", 
                 "Jun", "Lmna", "Myadm", "Nppc", "Phlda1", "Rhob", "Slc10a6", "Stat3", "Tra2b", "Zfand5", 
                 "Kcne4", "Atf4", "Cebpb", "Ddx5", "Egr2", "Fosl2", "Hspa1a", "Ier2", "Junb", "Maff", "Myc", 
                 "Pnp", "Rhoh", "Slc38a2", "Tagln2", "Trib1", "Zfp36", "Bag3", "Cebpd", "Des", "Eif1", 
                 "Gadd45a", "Hspa1b", "Ier3", "Jund", "Mafk", "Myd88", "Odc1", "Pnrc1", "Ripk1", "Slc41a1", 
                 "Tiparp", "Tubb4b", "Zfp36l1", "Bhlhe40", "Cebpg", "Dnaja1", "Eif5", "Gcc1", "Hspa5", "Ifrd1", 
                 "Klf2", "Mcl1", "Nckap5l", "Osgin1", "Ppp1cc", "Sat1", "Socs3", "Tnfaip3", "Tubb6", "Zfp36l2", 
                 "Brd2", "Csrnp1", "Dnajb1", "Gem", "Hspa8", "Il6", "Klf4", "Midn", "Ncoa7", "Oxnad1", 
                 "Ppp1r15a", "Sbno2", "Sqstm1", "Tnfaip6", "Ubc", "Zyx")

mapIt <- function(mouseids, horg, morg, orth) {
  mouseg <- mapIds(morg, mouseids, "ENTREZID", "SYMBOL")
  mapped <- select(orth, keys = mouseg, columns = "Homo_sapiens", keytype = "Mus_musculus")
  names(mapped) <- c("Mus_egid", "Homo_egid")
  husymb <- select(horg, keys = as.character(mapped$Homo_egid), columns = "SYMBOL", keytype = "ENTREZID")
  return(data.frame(Mus_symbol = mouseids,
                    Mus_egid = mouseg,
                    Homo_egid = mapped$Homo_egid,
                    Homo_symbol = husymb$SYMBOL))
}

human_mapped <- mapIt(mouse_genes, org.Hs.eg.db, org.Mm.eg.db, Orthology.eg.db)
human_mapped$Homo_symbol[human_mapped$Mus_symbol == "Cxcl1"] <- "CXCL1"
human_mapped$Homo_symbol[human_mapped$Mus_symbol == "Mir22hg"] <- "MIR22HG"
human_mapped$Homo_symbol[human_mapped$Mus_symbol == "Cyr61"] <- "CCN1"
human_mapped$Homo_symbol[human_mapped$Mus_symbol == "Fam132b"] <- "ERFE"
human_mapped$Homo_symbol[human_mapped$Mus_symbol == "Mt1"] <- "MT1A"
human_mapped$Homo_symbol[human_mapped$Mus_symbol == "Mt1"] <- "MT1A"
human_mapped$Homo_symbol[human_mapped$Mus_symbol == "Ccrn4l"] <- "NOCT"
human_mapped$Homo_symbol[human_mapped$Mus_symbol == "Hspa1b"] <- "HSPA1B"

DIG <-  human_mapped$Homo_symbol
DIG_list <- list(diss_genes = DIG)

# TCR genes ----
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
all_genes <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol", "description"),
  mart = mart
)
tcr_genes <- all_genes[grep("^TRAV|^TRBV|^TRGV|^TRDV", all_genes$hgnc_symbol), ]
tcr_genes <- tcr_genes$hgnc_symbol

# Ribosomal Proteins genes ----
rp_genes <- all_genes[grep("^RP([0-9]+-|[LS])", all_genes$hgnc_symbol), ]
rp_genes <- rp_genes$hgnc_symbol

#### ---- Representative Pre-processing and integration steps for cells of one macro-group in scRNAseq samples from one STS histology ---- ####
# load histology-specific integrated object with all cells and all samples 
UPS <- readRDS(".../UPS_object_integrated.rds")

#subset for one macro-group of interest - in this case T-Lymphocytes 
UPS_lymphoid <- subset(UPS, subset = classification == "Lymphoid Cells (TME)")

#remove cells that express HB genes (possible contaminants / RBCs)
hemoglbin_genes <- c("HBB", "HBA1", "HBA2", "HBBP1", "HBAP1", "HBD", 
                     "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ", "HBZP1")
hb_expression <- FetchData(UPS_lymphoid, 
                           vars = hemoglbin_genes, 
                           slot = "counts")
contaminated_cells <- rownames(hb_expression)[apply(hb_expression, 1, 
                                                    function(x) any(x != 0))]

UPS_lymphoid <- UPS_lymphoid[, !(colnames(UPS_lymphoid) %in% contaminated_cells)]

#Isolate B-cells
b_cells_genes <- c("CD19", "MS4A1")
b_cell_expression <- FetchData(UPS_lymphoid, 
                               vars = b_cells_genes, 
                               slot = "counts")
b_cells_tags <- rownames(b_cell_expression)[apply(b_cell_expression, 1, function(x) any(x != 0))]
#Subset for b-cells only 
UPS_b_cells <- UPS_lymphoid[, colnames(UPS_lymphoid) %in% b_cells_tags]

#Remove B-cells from T-cells object and analyze T-cells
UPS_lymphoid <- UPS_lymphoid[, !(colnames(UPS_lymphoid) %in% b_cells_tags)]

#Split by sample
UPS_lymphoid <- SplitObject(UPS_lymphoid, split.by = "sample")
# Calculate cell cycle score
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
UPS_lymphoid <- lapply(X = UPS_lymphoid, 
                       FUN = CellCycleScoring, 
                       s.features = s.genes, 
                       g2m.features = g2m.genes, 
                       set.ident = FALSE)
#Add DIG signature using AddModuleScore 
UPS_lymphoid <- lapply(X = UPS_lymphoid, 
                       FUN = AddModuleScore,
                       name = "DIG_score",
                       features = DIG_list)
#Run SCTransform and regress out cell cycle scores and mito genes pct
UPS_lymphoid <- lapply(X = UPS_lymphoid, 
                       FUN = SCTransform,
                       method = "glmGamPoi", 
                       vst.flavor = "v2",
                       vars.to.regress = c("percent.mt", 
                                           "S.Score", 
                                           "G2M.Score", 
                                           "DIG_score1"),
                       return.only.var.genes = FALSE)

var.features <- SelectIntegrationFeatures(object.list = UPS_lymphoid, nfeatures = 3000)
UPS_lymphoid <- merge(x = UPS_lymphoid[[1]], y = UPS_lymphoid[2:length(UPS_lymphoid)], merge.data=TRUE)

# Define blacklist gene list and remove them from from the var.features 
blacklist <- unique(c(DIG, tcr_genes, , "MALAT1", 
                      rp_genes))

var.features <- setdiff(var.features, blacklist)

VariableFeatures(UPS_lymphoid) <- var.features

UPS_lymphoid <- RunPCA(UPS_lymphoid, features = var.features,  verbose = TRUE)

UPS_lymphoid <- RunHarmony(UPS_lymphoid, assay.use="SCT", dims.use = 1:15,
                           group.by.vars = "sample")

UPS_lymphoid <- RunUMAP(UPS_lymphoid, reduction = "harmony", dims = 1:15)

UPS_lymphoid <- FindNeighbors(UPS_lymphoid, 
                              reduction = "harmony", 
                              k.param = 20,
                              dims = 1:15) %>% 
  FindClusters(resolution = 1.0) %>% 
  FindClusters(resolution = 1.1) %>% 
  FindClusters(resolution = 1.2) %>% 
  FindClusters(resolution = 1.3) %>% 
  FindClusters(resolution = 1.4) %>% 
  FindClusters(resolution = 1.5) %>% 
  FindClusters(resolution = 1.6) %>% 
  FindClusters(resolution = 1.7) %>% 
  FindClusters(resolution = 1.8) %>% 
  FindClusters(resolution = 1.9) %>% 
  FindClusters(resolution = 2.0) %>% 
  FindClusters(resolution = 0.9) %>% 
  FindClusters(resolution = 0.8) %>% 
  FindClusters(resolution = 0.7) %>% 
  FindClusters(resolution = 0.6) %>% 
  FindClusters(resolution = 0.5) %>% 
  FindClusters(resolution = 0.4)   

# Run clustree to empirically decide on the optimal clustering resolution
clustree(UPS_lymphoid, prefix = "SCT_snn_res.")

# Based on clustree res set to 1.3 
UPS_lymphoid <- FindNeighbors(UPS_lymphoid, 
                              reduction = "harmony", 
                              k.param = 20,
                              dims = 1:15) %>% 
  FindClusters(resolution = 1.3)

# Run Seurat's FindAllMarkers
ups_lymph_markers_sct_harmony <- FindAllMarkers(UPS_lymphoid, 
                                                min.pct = 0.25, 
                                                slot = "counts",
                                                assay = "RNA",
                                                only.pos = T)

#### ---- Representative code used to run limma-voom differential gene expression for marker gene selection ----
# Run limma-voom for DGE on macro-group specific Seurat objects 
# "myelo_all_sarc" is the myeloid cell object including all samples 
# "cell_id" contains cluster assignment information 

myelo_all_sarc$cell_id <- myelo_all_sarc$cell_id
myelo_all_sarc$cell_id <- gsub("_", "", myelo_all_sarc$cell_id)
myelo_all_sarc$cell_id <- gsub("[_+\\-]", "", myelo_all_sarc$cell_id)

myelo_all_sarc@meta.data$sampl_cell <- paste(myelo_all_sarc$sample, 
                                             myelo_all_sarc$cell_id, sep = "_")

myeloid_pseudobulk <- AggregateExpression(myelo_all_sarc, slot = "counts", 
                                          assays = "SCT", group.by = "sampl_cell")

myeloid_pseudobulk <- myeloid_pseudobulk$SCT

#Create DGEList object
dge <- DGEList(counts=myeloid_pseudobulk)
dge <- calcNormFactors(dge)

#filter for low expressed genes
cutoff <- 10
drop <- which(apply(cpm(dge), 1, max) < cutoff)
dge_filtered <- dge[-drop,] 
dim(dge_filtered)

matrix_metadata <- data.frame(sample=names(as.data.frame(dge_filtered)))
rownames(matrix_metadata) <- matrix_metadata$sample

matrix_metadata <- matrix_metadata %>%
  separate(col = sample, into = c("sample", "cell_id"), sep = "_", fill = "right")

# Create the design matrix
design <- model.matrix(~0 + matrix_metadata$cell_id)
colnames(design) <- levels(as.factor(matrix_metadata$cell_id))

# Correlation (duplicateCorrelation) by sample 
v <- voom(dge_filtered, design, plot = TRUE)
corfit <- duplicateCorrelation(v, design, block = matrix_metadata$sample)

# Re-run voom with correlation 
v <- voom(dge_filtered, design, plot = T, block = matrix_metadata$sample, 
          correlation = corfit$consensus)

fit <- lmFit(v, design)
fit <- eBayes(fit)

clusters <- colnames(coef(fit))

contr.list <- lapply(seq_along(clusters), function(i) {
  others <- clusters[-i]
  this.contrast <- paste0(
    clusters[i],
    " - (",
    paste(others, collapse = " + "),
    ")/",
    length(others)
  )
  names(this.contrast) <- paste0(clusters[i], ".vs.all")
  return(this.contrast)
})

contr.list <- unlist(contr.list)

contr.matrix <- makeContrasts(
  contrasts = contr.list,
  levels = design
)

fit2 <- contrasts.fit(fit, contr.matrix)
fit2 <- eBayes(fit2)

# Initialize an empty data frame for storing results
results_limma_myeloid <- data.frame(
  logFC = numeric(),
  AveExpr = numeric(),
  t = numeric(),
  P.Value = numeric(),
  adj.P.Val = numeric(),
  B = numeric(),
  cell_type = character(),
  gene = character(),
  logp = numeric(),
  stringsAsFactors = FALSE
)

# Loop over each contrast in contr.matrix
for (cn in colnames(contr.matrix)) {
  res <- topTable(fit2, coef = cn, number = Inf, sort.by = "P")
  res$gene <- rownames(res)
  rownames(res) <- NULL
  cell_type_name <-sub(" -.*", "", cn)
  res$cell_type <- cell_type_name
  res$logp <- -log10(res$adj.P.Val)
  results_limma_myeloid <- rbind(results_limma_myeloid, res)
}

