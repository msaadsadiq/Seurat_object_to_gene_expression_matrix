

####### ONE TIME FUNCTIONS #######

####### Make sure you are System Administrator on the PC you are-
####### performing the analysis on.
####### Make a new R-Project file before continuing :D

####### INITIALLIZING R
####### Setting library to your R folder in C. This is done because on-
####### installing gglot2 somehow R didnt know th Lib location

####### yes the installation still says (as ‘lib’ is unspecified) but dont-
####### worry, its working fine as long as its installing in your-
####### defined .libPath

.libPaths("C:/Users/mxq52/Documents/R/win-library/3.4")


####### Install Seurat v3.0 Prerelease Available here: 
####### https://satijalab.org/seurat/install.html#prerelease
####### Multiple errors appeared in the installation cycle, you must first- 
####### install ggplot2, cowplot and ggridges

install.packages('ggplot2')
install.packages("cowplot")
install.packages("ggridges")
install.packages('devtools')
install.packages('Seurat')

####### START FROM HERE EVERYTIME YOU OPEN THE PROJECT #######
####### Install these packages in the correct sequence, some will load from-
####### ggplot2 thats okay.

library(ggplot2)
library(cowplot)
library(Matrix)
library(ggridges)
library(Seurat)
library(dplyr)

####### MAKING R LOAD DATA #######
####### Read data onto R this will take a few minutes based on your system
alk3_1.data <- Read10X(data.dir = "C:/Users/mss255/Downloads/Fahd/01_analysis/ALK3-Tube1/filtered_gene_bc_matrices/GRCh38")



### Extract raw matrix ######################

alk3_1.df <- matrix(0, ncol=910, nrow=33694)
alk3_1.df <- as.data.frame(alk3_1.df)
colnames(alk3_1.df) <- alk3_1.data@Dimnames[[2]]
rownames(alk3_1.df) <- alk3_1.data@Dimnames[[1]]

index_final <- matrix(0, ncol=910, nrow=33694)
index <- alk3_1.data@i
index_j <- alk3_1.data@j
index<- as.data.frame(index)
index_j<- as.data.frame(index_j)

index_j$values <- c(NA)
index_j$index_j <- index_j$index_j+1

for(i in 1:910){
  index_j$values[index_j$index_j==i] <- index$index[index_j$index_j==i]
}

  for(i in 1:910){
    alk3_1.df[index_j$values[index_j$index_j==i],i] <- alk3_1.data@x[index_j$index_j==i]  
  }
  
  alk3_1.df_t <- t(alk3_1.df)
  alk3_1.df_t <- as.data.frame(alk3_1.df_t)
  
  alk3_1.df_t$human <- c("none")
  for(i in 1:910){
    alk3_1.df_t$human[i] <- paste("normal.human.1.",i,sep="")
  }
  
  alk3_1.df_t$barcode <- rownames(alk3_1.df_t)
  alk3_1.df_t$assigned_cluster <- 0 
  
  alk3_1.df_t <- alk3_1.df_t[,c(33695:33697, order(colnames(alk3_1.df_t[1:33694])))]
    


write.csv(file="rawMatrix.csv", alk3_1.df_t)


####### Examine memory savings between regular and sparse matrices
dense.size <- object.size(x = as.matrix(x = alk3_1.data)
                          )

dense.size

sparse.size <- object.size(x = alk3_1.data
                           )

sparse.size

dense.size/sparse.size

####### DATA INITIALIZATION #######
####### Initializing the ALK3 data set (naming based on the n, this value 
####### changes, e.g. alk3_1, alk3_2, alk3_3 etc)
####### We are selecting >= 1 cell since that represents 0.1% of data, having 
####### less than 200 detected genes.
alk3_1 <- CreateSeuratObject(raw.data = alk3_1.data, 
                             min.cells = 1, 
                             min.genes = 200,
                             project = "10X_alk3_1"
                             )

####### QC by mitochondrial genes is a common QC method in scRNAseq analysis
####### here we perform mito percentage based gene normalization
mito.genes <- grep(pattern = "^MT-", x = rownames(x = alk3_1@data), value = TRUE
                   )

percent.mito <- Matrix::colSums(alk3_1@raw.data[mito.genes, ])/Matrix::colSums(alk3_1@raw.data
                                                                               )

####### AddMetaData adds columns to object@meta.data, and is a great place to
####### stash QC stats
alk3_1 <- AddMetaData(object = alk3_1, 
                      metadata = percent.mito, 
                      col.name = "percent.mito"
                      )

VlnPlot(object = alk3_1, features.plot = c("nGene",
                                           "nUMI",
                                           "percent.mito"),
        nCol = 3
        )


####### Since there is a rare subset of cells with an outlier level of high 
####### mitochondrial percentage and also low UMI content, we filter- 
####### these as well, the first step is visualization
par(mfrow = c(1, 2)
    )

GenePlot(object = alk3_1, gene1 = "nUMI",
         gene2 = "percent.mito"
         )

GenePlot(object = alk3_1, gene1 = "nUMI",
         gene2 = "nGene"
         )

####### now that we can see the percentage of mitochondiral and gene outliers-
####### we select to remove CELLS >0.05% (total RNA) mitochondrial genes and-
####### between 2200 to 7000 unique gene counts/cell
alk3_1 <- FilterCells(object = alk3_1, 
                      subset.names = c("nGene",
                                       "percent.mito"),
                      low.thresholds = c(2200, -Inf),
                      high.thresholds = c(7000, 0.05)
                      )


####### DATA NORMALIZATION #######
####### Scales by x10,000 by default and then log norms it
alk3_1 <- NormalizeData(object = alk3_1, 
                        normalization.method = "LogNormalize", 
                        scale.factor = 10000
                        )


####### HIGHLY VARIABLE GENES #######
alk3_1 <- FindVariableGenes(object = alk3_1, mean.function = ExpMean, 
                            dispersion.function = LogVMR, 
                            x.low.cutoff = 0.0125, 
                            x.high.cutoff = 3, 
                            y.cutoff = 0.5
                            )

length(x = alk3_1@var.genes)


####### DATA SCALLING #######
alk3_1 <- ScaleData(object = alk3_1, 
                    vars.to.regress = c("nUMI", 
                                       "percent.mito")
                    )


####### LINEAR DIMENSIONAL REDUCTION #######
alk3_1 <- RunPCA(object = alk3_1, 
                 pc.genes = alk3_1@var.genes, 
                 do.print = TRUE, 
                 pcs.print = 1:5, 
                 genes.print = 5
                 )


####### DATA VISUALIZATION #######
####### data examination
####### Multi-PCA visualization
PrintPCA(object = alk3_1, pcs.print = 1:5, 
         genes.print = 5, 
         use.full = FALSE
         )


VizPCA(object = alk3_1, 
       pcs.use = 1:4
       )


PCAPlot(object = alk3_1, 
        dim.1 = 1, 
        dim.2 = 2, 
        use.full = T
        )


####### ProjectPCA scores each gene in the dataset (including genes not 
####### included in the PCA) based on their correlation with the 
####### calculated components.It can be used to identify markers
####### that are strongly correlated with cellular heterogeneity, 
####### but may not have passed through variable gene selection.
####### The results of the projected PCA can be explored by setting 
####### use.full=T in the functions above
alk3_1 <- ProjectPCA(object = alk3_1, 
                     do.print = FALSE
                     )

PCHeatmap(object = alk3_1, 
          pc.use = 1, 
          cells.use = 500, 
          do.balanced = TRUE, 
          label.columns = FALSE
          )


####### splitting PCA into 12 clusters
PCHeatmap(object = alk3_1, 
          pc.use = 1:10, 
          cells.use = 500, 
          do.balanced = TRUE, 
          label.columns = FALSE, 
          use.full = FALSE
          )

####### STASTICALLY SIGNIFICANT PCA #######
####### anything below the jackstrw (low pvalues) provides an idea on PC number
alk3_1 <- JackStraw(object = alk3_1, num.replicate = 100, 
                    display.progress = TRUE
                    )

JackStrawPlot(object = alk3_1, 
              PCs = 1:12
              )

####### alternative method using PCElbow (computational conservaion)
PCElbowPlot(object = alk3_1
            )


####### CLUSTERING CELLS #######
####### save.SNN = T saves the SNN so that the clustering algorithm can be rerun
####### force.recalc forces a rerun and resolution increases from 0.6-1.2
####### more resolution more clusters
alk3_1 <- FindClusters(object = alk3_1, 
                       reduction.type = "pca", 
                       dims.use = 1:10,
                       resolution = 0.6, 
                       print.output = 0, 
                       save.SNN = TRUE,
                       force.recalc = TRUE
                       )

PrintFindClustersParams(object = alk3_1
                        )


####### NON LINEAR DIMENSIONALITY REDUCTION tSNE #######
alk3_1 <- RunTSNE(object = alk3_1, 
                  dims.use = 1:10, 
                  do.fast = TRUE
                  )

TSNEPlot(object = alk3_1
         )

####### saving the tSNE plot
saveRDS(alk3_1, file = "C:/Users/mxq52/Box/scRNAseq Analysis HIRN/HP2300/scRNAseq R project/ALK3bright1/alk3_1_tSNE.rds")

###### REALM OF GENES #######
###### DIFFERENTIALLY EXPRESSED GENE ANALYSIS #######
###### classifying all markers of cluster0
cluster0.markers <- FindMarkers(object = alk3_1, 
                                ident.1 = 0, 
                                min.pct = 0,
                                thresh.use = 0
                                )

print(x = head(x = cluster0.markers, 
               n = 5)
      )

###### classifying all markers of cluster1
cluster1.markers <- FindMarkers(object = alk3_1, 
                                ident.1 = 1, 
                                min.pct = 0,
                                thresh.use = 0
                                )

print(x = head(x = cluster1.markers, 
               n = 5)
      )

###### classifying all markers of cluster2
cluster2.markers <- FindMarkers(object = alk3_1, 
                                ident.1 = 2, 
                                min.pct = 0,
                                thresh.use = 0
                                )

print(x = head(x = cluster2.markers, 
               n = 5)
      )

###### classifying all markers of cluster3
cluster3.markers <- FindMarkers(object = alk3_1, 
                                ident.1 = 3, 
                                min.pct = 0,
                                thresh.use = 0
                                )

print(x = head(x = cluster3.markers, 
               n = 5)
      )

###### classifying all markers of cluster4
cluster4.markers <- FindMarkers(object = alk3_1, 
                                ident.1 = 4, 
                                min.pct = 0,
                                thresh.use = 0
                                )

print(x = head(x = cluster4.markers, 
               n = 5)
      )

###### classifying all markers of cluster5
cluster5.markers <- FindMarkers(object = alk3_1, 
                                ident.1 = 5, 
                                min.pct = 0,
                                thresh.use = 0
                                )

print(x = head(x = cluster5.markers, 
               n = 5)
      )

####### an example of a VERY stringent test ROC classification

#cluster1.markers <- FindMarkers(object = alk3_1, 
#                               ident.1 = 0, 
#                                thresh.use = 0.25, 
#                                test.use = "roc", 
#                                only.pos = TRUE
#                                )


###### find markers for every cluster compared to all remaining cells, report
###### only the positive ones, the difference to this and that above is this
###### is for all clusters and is fairly faster for older versions of seurat
###### this is called find_all_markers

# Define an order of cluster identities this works
my_levels <- c(2,3,0,1,4)

# Relevel object@ident
alk3_1@ident <- factor(x = alk3_1@ident, levels = my_levels)

alk3_1.markers <- FindAllMarkers(object = alk3_1, 
                                 only.pos = TRUE, 
                                 min.pct = 0,
                                 thresh.use = 0
                                 )

alk3_1.markers %>% group_by(cluster) %>% top_n(2, avg_logFC
                                               )



####### Plotting data using violin plots
VlnPlot(object = alk3_1, features.plot = c("KIF1A",
                                           "MAP1A",
                                           "SFRP5",
                                           "FXYD2",
                                           "FGR",
                                           "GPX2"), 
        y.max = 4
        )


####### Using raw UMI counts
VlnPlot(object = alk3_1, features.plot = c("KIF1A",
                                           "MAP1A",
                                           "SFRP5",
                                           "FXYD2",
                                           "FGR",
                                           "GPX2"),
        use.raw = TRUE, 
        y.log = TRUE
        )


####### Using raw UMI counts for Acinar genes
VlnPlot(object = alk3_1, features.plot = c ("CTRB1", 
                                            "CELA3A", 
                                            "CELA3B", 
                                            "CTRB2", 
                                            "PLA2G1B", 
                                            "PRSS2", 
                                            "SPINK1", 
                                            "CLPS", 
                                            "CPA1", 
                                            "PRSS1", 
                                            "CPA2", 
                                            "REG1A",
                                            "PNLIP",
                                            "SYCN",
                                            "PNLIPRP1",
                                            "CTRC",
                                            "KLK1",
                                            "CELA2A",
                                            "CPB1"),
        use.raw = TRUE, 
        y.log = TRUE
)

####### Using raw UMI counts for ductal genes
VlnPlot(object = alk3_1, features.plot = c ("KRT19", 
                                            "TACSTD2", 
                                            "ANXA2", 
                                            "S100A10", 
                                            "S100A11", 
                                            "KRT17", 
                                            "KRT18", 
                                            "KRT7", 
                                            "S100A16", 
                                            "S100A14", 
                                            "TPM1" 
                                            ),
        use.raw = TRUE, 
        y.log = TRUE
)

####### Using raw UMI counts for ductal genes
VlnPlot(object = alk3_1, features.plot = c ("ONECUT2",
                                            "LITAF",
                                            "SOX4",
                                            "DAB2",
                                            "CREB5",
                                            "HLA-DQB1",
                                            "WWTR1",
                                            "PPARGC1A",
                                            "PKHD1",
                                            "NFIB"
                                            ),
use.raw = TRUE, 
y.log = TRUE
)

####### make a tsne of multiple genes in one panel

FeaturePlot(object = alk3_1, features.plot = c("KRT19", 
                                               "COL1A1", 
                                               "BMPR1A", 
                                               "PDX1", 
                                               "SOX9", 
                                               "CA2", 
                                               "KIF1A", 
                                               "MAP1A", 
                                               "SFRP5", 
                                               "FXYD2", 
                                               "FGR", 
                                               "GPX2"), 
            cols.use = c("light grey", 
                         "red"), 
            reduction.use = "tsne"
            )

####### make a tsne of multiple genes in one panel

FeaturePlot(object = alk3_1, features.plot = c("KRT19", 
                                               "COL1A1", 
                                               "BMPR1A", 
                                               "PDX1", 
                                               "SOX9", 
                                               "CA2", 
                                               "KIF1A", 
                                               "MAP1A", 
                                               "SFRP5", 
                                               "FXYD2", 
                                               "FGR", 
                                               "GPX2"), 
            cols.use = c("light grey", 
                         "red"), 
            reduction.use = "tsne"
)

####### make a tsne of multiple genes in one panel

FeaturePlot(object = alk3_1, features.plot = c("CTRB1", 
                                               "CELA3A", 
                                               "CELA3B", 
                                               "CTRB2", 
                                               "PLA2G1B", 
                                               "PRSS2", 
                                               "SPINK1", 
                                               "CLPS", 
                                               "CPA1", 
                                               "PRSS1", 
                                               "CPA2", 
                                               "REG1A",
                                               "PNLIP",
                                               "SYCN",
                                               "PNLIPRP1",
                                               "CTRC",
                                               "KLK1",
                                               "CELA2A",
                                               "CPB1"
                                               ), 
            cols.use = c("light grey", 
                         "red"), 
            reduction.use = "tsne"
)

####### setting slim.col.label to TRUE will print just the cluster IDS instead of
####### every cell name
####### Heatmap of 6 clusters, naming is from 0-5 so its a good idea to look at
####### the cluster from earlier to be sure regarding cluster assignment. 
####### Clustering is still unbiased, you are just re-classifying them now
####### Color Pallette https://htmlcolorcodes.com/
top10 <- alk3_1.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)

DoHeatmap(object = alk3_1, 
          genes.use = top10$gene, 
          slim.col.label = TRUE, 
          remove.key = TRUE,
          cells.use = WhichCells(alk3_1,c (0,1,2,3,4,5)),
          col.low = "#004DFF",
          col.mid = "#FFFECE",
          col.high = "#CF0000",
          rotate.key = FALSE,
          title = NULL,
          cex.col = 10,
          cex.row = 10,
          group.label.loc = "bottom",
          group.label.rot = FALSE,
          group.cex = 15,
          group.spacing = 0.15,
          assay.type = "RNA",
          do.plot = TRUE
)

####### Assigning cell type identity to clusters


current.cluster.ids <- c(0, 
                         1, 
                         2, 
                         3, 
                         4, 
                         5
                         )

new.cluster.ids <- c("GPX2+ Ductal", 
                     "KIF1A+ Ductal", 
                     "COL1A1+ Mesenchymal", 
                     "SFRP5+ Ductal", 
                     "CA2+ Ductal", 
                     "FGR+ Ductal?"
                     )

alk3_1@ident <- plyr::mapvalues(x = alk3_1@ident, 
                                from = current.cluster.ids, 
                                to = new.cluster.ids)

TSNEPlot(object = alk3_1, 
         do.label = TRUE, 
         pt.size = 2) 

####### Further subdivisions within cell types
####### First lets stash our identities for later

alk3_1 <- StashIdent(object = alk3_1, 
                     save.name = "ClusterNames_0.6")

alk3_1 <- FindClusters(object = alk3_1, reduction.type = "pca", dims.use = 1:10, 
                     resolution = 1.2, print.output = FALSE)

####### How to plot two tSNE plots side by side, and how to color
####### points based on different criteria


plot1 <- TSNEPlot(object = alk3_1, 
                  do.return = TRUE, 
                  no.legend = TRUE,
                  pt.size = 2.5,
                  do.label = TRUE)

plot2 <- TSNEPlot(object = alk3_1, 
                  do.return = TRUE, 
                  group.by = "ClusterNames_0.6", 
                  no.legend = TRUE,
                  pt.size = 2.5,
                  do.label = TRUE)

plot_grid(plot1, plot2)

####### Find discriminating markers
ductcell.markers <- FindMarkers(object = alk3_1, ident.1 = 0, ident.2 = 1)

FeaturePlot(object = alk3_1, features.plot = c("OLFM4", "REG1B"), cols.use = c("yellow", 
                                                                             "red"))
alk3_1 <- SetAllIdent(object = alk3_1, id = "ClusterNames_0.6")
saveRDS(alk3_1, file = "C:/Users/mxq52/Box/scRNAseq Analysis HIRN/HP2300/scRNAseq R project/ALK3bright1ALK3bright1final.rds")









