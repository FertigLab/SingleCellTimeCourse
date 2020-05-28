###################################################
# code chunk: init
###################################################
library(cellrangerRkit) # to open cellranger files
library(devtools)
library(monocle)
library(reticulate)
library(Matrix)
library(VGAM) # for the `expressionFamily` function when creating the `CellDataSet` object
library(clues)
library(dplyr)
library(mcclust)
library(flexclust)
library(DelayedArray)
import("louvain") # import louvain before umap, otherwise it will not find the module
import("umap")
library(ggplot2)

###################################################
# code chunk: load and preparing data
###################################################

# loading CellRanger data
load("CellRanger_LTK.RData")

# SAMPLE IDS:
# LTK1C - SCC1 CTX REP1
# LTK1CR - SCC1 CTX REP2
# LTK6C - SCC6 CTX REP1
# LTK6CR - SCC6CTX REP2
# LTK25C - SCC25 CTX REP1
# LTK25CR - SCC25 CTX REP2
# LTK1P - SCC1 PBS REP1
# LTK1PR - SCC1 PBS REP2
# LTK6P - SCC6 PBS REP1
# LTK6PR - SCC6 PBS REP2
# LTK25P - SCC25 PBS REP1
# LTK25PR - SCC25 PBS REP2

### Formating and creating matrices to merge

## SCC1-CTX1
SCC1.CTX1.mat <- as.matrix(exprs(LTK1C))
# check matrix
SCC1.CTX1.mat[1:5,1:5]
# check dims
dim(SCC1.CTX1.mat)
# read in gene names 
gene.symbols <- fData(LTK1C)$symbol
# set gene short names to rownames of matrix
rownames(SCC1.CTX1.mat) <- gene.symbols
# check new rownames
SCC1.CTX1.mat[1:5,1:5]
# read in sample names
SCC1.CTX1.barcodes <- pData(LTK1C)$barcode
# set barcodes to colnames of matrix
colnames(SCC1.CTX1.mat) <- SCC1.CTX1.barcodes
# check new colnames
SCC1.CTX1.mat[1:5,1:5]
# add unique identifier to colnames
colnames(SCC1.CTX1.mat) <- paste(colnames(SCC1.CTX1.mat), "SCC1", sep = "_")
colnames(SCC1.CTX1.mat) <- paste(colnames(SCC1.CTX1.mat), "CTX", sep = "_")
colnames(SCC1.CTX1.mat) <- paste(colnames(SCC1.CTX1.mat), "R1", sep = "_")
# pData for treatment
SCC1.CTX1.pdat <- data.frame("samples" = colnames(SCC1.CTX1.mat), "cell" = "SCC1", "treatment" = "CTX", "replicate" = "R1")

### Repeat the steps above for all samples

### Join data to create CDS

# join matrices
exprsMat <- cbind(SCC1.CTX1.mat,SCC1.CTX2.mat,SCC1.PBS1.mat,SCC1.PBS2.mat,SCC6.CTX1.mat,SCC6.CTX2.mat,SCC6.PBS1.mat,SCC6.PBS2.mat,SCC25.CTX1.mat,SCC25.CTX2.mat,SCC25.PBS1.mat,SCC25.PBS2.mat)
save(exprsMat, file="exprsMat.Rda")

exprsSCC1Mat <- cbind(SCC1.CTX1.mat,SCC1.CTX2.mat,SCC1.PBS1.mat,SCC1.PBS2.mat)
save(exprsSCC1Mat, file="exprsSCC1Mat.Rda")

exprsSCC6Mat <- cbind(SCC6.CTX1.mat,SCC6.CTX2.mat,SCC6.PBS1.mat,SCC6.PBS2.mat)
save(exprsSCC6Mat, file="exprsSCC6Mat.Rda")

exprsSCC25Mat <- cbind(SCC25.CTX1.mat,SCC25.CTX2.mat,SCC25.PBS1.mat,SCC25.PBS2.mat)
save(exprsSCC25Mat, file="exprsSCC25Mat.Rda")

# check that dim is the same
dim(exprsMat)
dim(exprsSCC1Mat)
dim(exprsSCC6Mat)
dim(exprsSCC25Mat)

## make pData dictionary
# join pData
pdat <- as.data.frame(rbind(SCC1.CTX1.pdat,SCC1.CTX2.pdat,SCC1.PBS1.pdat,SCC1.PBS2.pdat,SCC6.CTX1.pdat,SCC6.CTX2.pdat,SCC6.PBS1.pdat,SCC6.PBS2.pdat,SCC25.CTX1.pdat,SCC25.CTX2.pdat,SCC25.PBS1.pdat,SCC25.PBS2.pdat))
rownames(pdat) <- pdat$samples
save(pdat, file="phenoData.Rda")

pdatSCC1 <- as.data.frame(rbind(SCC1.CTX1.pdat,SCC1.CTX2.pdat,SCC1.PBS1.pdat,SCC1.PBS2.pdat))
rownames(pdatSCC1) <- pdatSCC1$samples
save(pdatSCC1, file="phenoDataSCC1.Rda")

pdatSCC6 <- as.data.frame(rbind(SCC6.CTX1.pdat,SCC6.CTX2.pdat,SCC6.PBS1.pdat,SCC6.PBS2.pdat))
rownames(pdatSCC6) <- pdatSCC6$samples
save(pdatSCC6, file="phenoDataSCC6.Rda")

pdatSCC25 <- as.data.frame(rbind(SCC25.CTX1.pdat,SCC25.CTX2.pdat,SCC25.PBS1.pdat,SCC25.PBS2.pdat))
rownames(pdatSCC25) <- pdatSCC25$samples
save(pdatSCC25, file="phenoDataSCC25.Rda")

# generate fData
fdat <- toupper(as.matrix(gene.symbols))
# set gene short names as rownames
rownames(fdat) <- gene.symbols
fdat <- as.data.frame(fdat)
fdat$feature <- fdat
fdat$feature <- "NA"
colnames(fdat) <- c("gene_short_name", "feature")
save(fdat, file="featureData.Rda")

###################################################
# code chunk: Monocle - creating CellDataSet and normalization
###################################################

##### now the 3 datasets can be loaded to create CellDataSet objects
load("exprsMat.Rda")
load("featureData.Rda")
load("phenoData.Rda")
colnames(exprsMat) <- rownames(pdat)
rownames(exprsMat) <- rownames(fdat)

### generate CDS
cds.stc <- newCellDataSet(exprsMat,
                          phenoData = new("AnnotatedDataFrame", data = pdat),
                          featureData = new("AnnotatedDataFrame", data = fdat),
                          lowerDetectionLimit = 1,
                          expressionFamily = negbinomial.size())
save(cds.stc, file = "ShortTimeCourse_cds.Rda")

### estimate size factors and dispersions (normalization)
DelayedArray:::set_verbose_block_processing(TRUE)
options(DelayedArray.block.size=1000e6)
cds.stc <- estimateSizeFactors(cds.stc)
cds.stc <- estimateDispersions(cds.stc)

### reduce dimensionality
disp.table <-  dispersionTable(cds.stc)
disp.table <-  disp.table %>% mutate(excess_disp = (dispersion_empirical - dispersion_fit) / dispersion_fit) %>% arrange(plyr::desc(excess_disp))
top.genes <- as.character(head(disp.table, 2500)$gene_id)
cds.stc <- setOrderingFilter(cds.stc, top.genes)
plot_ordering_genes(cds.stc)
plot_pc_variance_explained(cds.stc, return_all = FALSE)

### normalization
cds.stc <- preprocessCDS(cds.stc, method = "PCA", norm_method = "log", num_dim = 50, verbose = TRUE)
cds.stc <- reduceDimension(cds.stc, max_components = 2, num_dim = 5, reduction_method = "tSNE", metric = "correlation", min_dist = 0.75, n_neighbors = 50, verbose = TRUE)

###################################################
# code chunk: Monocle - cell clusters (Figure 1)
###################################################

### clustering by cell line
plot_cell_clusters(cds.stc, 
                   color_by = 'cell', 
                   cell_size = 0.3, 
                   show_group_id = T) + 
  theme(legend.text=element_text(size=10)) + 
  theme(legend.position="right")

### clustering by treatment
plot_cell_clusters(cds.stc, 
                   color_by = 'treatment', 
                   cell_size = 0.3, show_group_id = T) + 
  scale_color_manual(breaks = c("CTX", "PBS"), values=c("red", "black")) +
  theme(legend.text=element_text(size=10)) + 
  theme(legend.position="right") 

### TAFP2a and VIM gene expression per cell
TFAP2A_id <- row.names(subset(fData(cds.stc), gene_short_name == "TFAP2A"))
VIM_id <- row.names(subset(fData(cds.stc), gene_short_name == "VIM"))
cth.ap2.vim <- newCellTypeHierarchy()
cth.ap2.vim <- addCellType(cth.ap2.vim, "TFAP2A+/VIM+", classify_func = function(x)  { x[TFAP2A_id,] > 1 & x[VIM_id,] > 2 })
cth.ap2.vim <- addCellType(cth.ap2.vim, "TFAP2A-/VIM-", classify_func = function(x)  { x[TFAP2A_id,] < 1 & x[VIM_id,] < 2 })
cth.ap2.vim <- addCellType(cth.ap2.vim, "TFAP2A+/VIM-", classify_func = function(x)  { x[TFAP2A_id,] > 1 & x[VIM_id,] < 2 })
cth.ap2.vim <- addCellType(cth.ap2.vim, "TFAP2A-/VIM+", classify_func = function(x)  { x[TFAP2A_id,] < 1 & x[VIM_id,] > 2 })
cds.stc <- classifyCells(cds.stc, cth.ap2.vim)
table(pData(cds.stc)$CellType)

### annotating cells according to TFAP2A and VIM expression
cds.stc <- clusterCells(cds.stc, method = 'louvain', res = 5e-4, louvain_iter = 1, verbose = T)
plot_cell_clusters(cds.stc, 
                   color_by = 'CellType', 
                   cell_size = 0.3, 
                   show_group_id = T) + 
  scale_color_manual(breaks = c("TFAP2A+/VIM+", "TFAP2A+/VIM-", "TFAP2A-/VIM+", "TFAP2A-/VIM-", "Unknown"), values=c("blue", "orange", "green", "purple", "gray")) +
  theme(legend.text=element_text(size=10)) + 
  theme(legend.position="right")