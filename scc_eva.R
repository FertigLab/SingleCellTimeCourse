###################################################
### code chunk: init
###################################################
library(GSReg)
library(Rmagic)
library(ggplot2)

setwd('/Users/emily/OneDrive - Johns Hopkins University/projects/Lu')
###################################################
### code chunk: load data
###################################################
## cell line data
load("exprsSCC1Mat.Rda")
load("exprsSCC6Mat.Rda")
load("exprsSCC25Mat.Rda")
## pathways
load("hallmark_pathways.rda")
load("TFAP2A_EMT.Rda")

###################################################
### code chunk: magic imputation - run in command line R
###################################################
## remove genes with standard dev of 0
## transpose expression matrices
removezeros <- function(p) {
print(sum(apply(p,1,sd)==0,na.rm=T))
keepIndex = (apply(p,1,sd)!=0)
#sum(keepIndex);sum(!keepIndex)
q = p[keepIndex,]
print(dim(q))
r <- t(q)
return(r)
}

scc1 <- removezeros(exprsSCC1Mat)
scc6 <- removezeros(exprsSCC6Mat)
scc25 <- removezeros(exprsSCC25Mat)

## impute 
mscc1 <- magic(scc1)
mscc6 <- magic(scc6)
mscc25 <- magic(scc25)

save(mscc1, file = "magic_scc1.rda")
save(mscc6, file = "magic_scc6.rda")
save(mscc25, file = "magic_scc25.rda")

load("magic_scc1.rda")
load("magic_scc6.rda")
load("magic_scc25.rda")

## get binary vectors for each cell line 
split1 <- strsplit(colnames(exprsSCC1Mat),split='_', fixed=TRUE)
c1 <- sapply(split1, "[", 3)
split2 <- strsplit(colnames(exprsSCC6Mat),split='_', fixed=TRUE)
c6 <- sapply(split2, "[", 3)
split3 <- strsplit(colnames(exprsSCC25Mat),split='_', fixed=TRUE)
c25 <- sapply(split3, "[", 3)

###################################################
### code chunk: EVA 
###################################################

do_eva <- function(data, type, pathway, pheno) {
  ## transpose back for EVA
  edata <- t(data[[1]])
  exprsdata <- as.matrix(edata)
  rownames(exprsdata) <- rownames(edata)
  
  phenobin <- as.integer(as.factor(pheno))
  ## classic eva
  VarAnKendallV = GSReg.GeneSets.EVA(geneexpres=exprsdata, pathways=pathway, phenotypes=factor(phenobin)) 
  pvalustat = sapply(VarAnKendallV,function(x) x$pvalue);
  kendall <- lapply(VarAnKendallV, function (x) x[1:20])
  dd  <-  as.data.frame(matrix(unlist(kendall), nrow=length(unlist(kendall[1]))))
  rownames(dd) <- make.unique(names(kendall[[1]]))
  colnames(dd) <- names(kendall)
  dd <- as.data.frame(t(dd))
  dd$pathway_name <- rownames(dd)
  
  ## add meta information
  dd$metric <- "kendall_tau"
  dd$cell_line <- type
  dd$group1 <- "CTX"
  dd$group2 <- "PBS"
  return(dd)
}

a1 <- do_eva(data = mscc1, type = "SCC1", pathway = hallmark_pathways, pheno = c1)
a2 <- do_eva(data = mscc1, type = "SCC1", pathway = TFAP2A_EMT, pheno = c1)
print("comparison 1 done")

b1 <- do_eva(data = mscc6, type = "SCC6", pathway = hallmark_pathways, pheno = c6)
b2 <- do_eva(data = mscc6, type = "SCC6", pathway = TFAP2A_EMT, pheno = c6)
print("comparison 2 done")

c1 <- do_eva(data = mscc25, type = "SCC25", pathway = hallmark_pathways, pheno = c25)
c2 <- do_eva(data = mscc25, type = "SCC25", pathway = TFAP2A_EMT, pheno = c25)
print("comparison 3 done")

save(a1,a2,b1,b2,c1,c2, file = "SCC_EVA.rda")


###################################################
### code chunk: boxplots
###################################################
load("scc25_EVAhallmark.Rda") # a5
load("scc25_EVATFAP2A_EMT.Rda") #a6
load("scc6_EVAhallmark.Rda") #a3
load("scc6_EVATFAP2A_EMT.Rda") #a4
load("scc1_EVAhallmark.Rda") #a1
load("scc1_EVATFAP2A_EMT.Rda") #a2

box_plots <- function(data, path) {
## mult test correction
data$adj.pval <- p.adjust(data$pvalue, method = "BH")
print("number of significant pathways")
print(nrow(data[data$adj.pval < 0.05,]))
## split to join in long format
E1 <- data[c(1,21,24)]
E2 <- data[c(2,21,25)]
colnames(E1) <- c("E","pathway","group")
colnames(E2) <- c("E","pathway","group")
plotdat <- rbind(E1,E2)

title <- paste("cell line", unique(data$cell_line), path)
p <- ggplot(plotdat, aes(x=group, y=E, fill = group)) +
  geom_boxplot(outlier.shape = NA) + theme_bw() +
  theme(axis.text.x=element_text(size=rel(0.75), angle=45, hjust = 1),
        axis.ticks.x=element_blank()) + ggtitle(title) +
  scale_fill_manual(values=c("blue", "red")) +
  geom_point(aes(fill = group, alpha = 0.8), size = 1.5, shape = 21, position = position_jitterdodge()) + 
  ylab("EVA statistic") + xlab("treatment")
return(p)
}

#tiff("cibersort_scores_NK.tiff", width = 6000, height = 5000, units = "px", res = 800)
pdf("EVA_boxplots_SCC.pdf")
box_plots(data = a1, path = "hallmark")
box_plots(data = a2, path = "EMT")
box_plots(data = a3, path = "hallmark")
box_plots(data = a4, path = "EMT")
box_plots(data = a5, path = "hallmark")
box_plots(data = a6, path = "EMT")
dev.off()

###################################################
### code chunk: heatmaps 
###################################################
library(ComplexHeatmap)
library(reshape2)

heatmaps <- function(data, path) {
## split to join in long format
E1 <- data[c(1,21,24)]
E2 <- data[c(2,21,25)]
colnames(E1) <- c("E","pathway","group")
colnames(E2) <- c("E","pathway","group")
plotdat <- rbind(E1,E2)

mat <- acast(plotdat, pathway ~ group , value.var='E')
mats <- t(apply(mat,2,scale))
rownames(mats) <- colnames(mat)

## add heatmap annotation
ha_column = columnAnnotation(foo2 = as.factor(rownames(mats)))
m <- Heatmap(t(mats), col=inferno(5), name = "mat", 
        top_annotation=ha_column, 
        #left_annotation=row_ha,
        border = TRUE,
        #show_column_names = FALSE
)

return(m)
}

pdf("EVA_heatmaps_scaled.pdf",width=15,height=10, paper = "USr")
heatmaps(data = a1, path = "hallmark")
heatmaps(data = a2, path = "EMT")
heatmaps(data = a3, path = "hallmark")
heatmaps(data = a4, path = "EMT")
heatmaps(data = a5, path = "hallmark")
heatmaps(data = a6, path = "EMT")
dev.off()

