#Functions
library(reticulate)
use_virtualenv("myenv")

library(scran)
library(scater)
library(cowplot)
library(pheatmap)
library(DropletUtils)
library(igraph)
library(SoupX)
library(pbapply)
library(paletteer)

#helpful themes
figure_theme <- theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
                      axis.line = element_blank(),
                      axis.text.x  = element_text(color = "black", size = 7.5),
                      axis.text.y  = element_text(color = "black", size = 7.5),
                      axis.title.x  = element_text(color = "black", size = 10),
                      axis.title.y  = element_text(color = "black", size = 10)
)
figure_colors <- paletteer::scale_color_paletteer_d("ggsci", "default_igv")
figure_fills <- paletteer::scale_fill_paletteer_d("ggsci", "default_igv")

# Get density of points in 2 dimensions.
#' @param x A numeric vector.
#' @param y A numeric vector.
#' @param n Create a square n by n grid to compute density.
#' @return The density within each square.
get_density <- function(x, y, n = 100) {
  require(MASS)
  dens <- kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

# Plot the x-y density with ggplot2 raster
#' @param x A numeric vector.
#' @param y A numeric vector.
#' @param n Create a square n by n grid to compute density.
#' @param xlab x axis label.
#' @param ylab y axis label.
#' @param title Title for the plot.
#' @return A plot of x and y colored by density.
density_plot_raster <- function(x, y, xlab = NULL, ylab = NULL, title = NULL, n = 100, colorscale = viridis::viridis(100)){
  require(ggrastr)
  dat <- data.frame(xlab = x, ylab = y)
  dat$density <- get_density(dat$x, dat$y, n=n)
  ggplot(dat, aes(x = x, y = y, col = density)) + ggrastr::geom_point_rast(pch=19, cex=0.05) + 
    scale_color_gradientn(colours = colorscale)+ xlab(xlab) + ylab(ylab) + ggtitle(title) + theme_bw() + theme(axis.text = element_text(color = "black"))
}

#plot out a feature on a ggplot - uses the viridis color scheme
#' @param layout - a layout to use, be it PCA, tSNE, UMAP, graph layout
#' @param gene - the gene symbol to use (character)
#' @param adata - the adata object to use
#' @param xlab - the x axis label
#' @param ylab - the y axis label
feature_plot <- function(layout = "X_umap", gene, adata, xlab = "Dim1", ylab="Dim2", size = 0.05, scale = FALSE){
  layout <- adata$obsm$get(layout)
  set.seed(100)
  idx <- sample(1:nrow(layout))
    g_idx <- which(adata$var$Symbol == gene)
  exprs_in <- t(adata$X)[g_idx, idx ]
  layout <- layout[idx, ]
  #this function plots the highest expressing cells on the top
  ggplot(data.frame("x" = layout[, 1], "y" = layout[, 2], "Expression" = exprs_in),
         aes(x = x, y=y, col = Expression)) + ggrastr::geom_point_rast(pch=19, cex=size) + xlab(xlab) + ylab(ylab)+
    scale_color_gradientn(colours = c("blue", "grey", "red")) + ggtitle(gene) + theme_classic()
}

#plot variance
#' @param adata - anndata object with PCA computed
plot_pc_variance <- function(adata){
  df <- data.frame("variance" = adata$uns$get("pca")$variance_ratio, "PC" = 1:length(adata$uns$get("pca")$variance_ratio))
  ggplot(df, aes(x = PC, y = variance)) + geom_point(pch = 21, size = 1.5, fill = "grey50") + theme_bw() + theme(axis.text = element_text(color = "black"))
}

#plot out an attribute on a ggplot - uses the viridis color scheme
#' @param layout - a layout to use, be it PCA, tSNE, UMAP, graph layout
#' @param attribute - the attribute to use
#' @param adata - the sce object to reference
#' @param xlab - the x axis label
#' @param ylab - the y axis label
attribute_plot <- function(layout = "X_umap", attribute, adata, xlab = "Dim1", ylab="Dim2", size = 0.05, log10_transform = FALSE){
  layout <- adata$obsm$get(layout)
  if(log10_transform == TRUE){
    attr_in <- log10(adata$obs[, attribute])
  }else{
    attr_in <- adata$obs[, attribute]
  }
  #this function plots the highest expressing cells on the top
  ggplot(data.frame("x" = layout[, 1], "y" = layout[, 2], "attr_in" = attr_in)[order(attr_in, decreasing = FALSE), ],
         aes(x = x, y=y, col = attr_in)) + ggrastr::geom_point_rast(pch=19, cex=size) + xlab(xlab) + ylab(ylab)+
    scale_color_viridis_c(1000) + ggtitle(attribute) + theme_classic() + theme(legend.title = element_blank())
}


#feature violin plots
#' @param adata - adaa object to use
#' @param gene - gene to plot
#' @param group - the grouping to use in adata$obs
feature_violin_plot <- function(adata, gene, group = "leiden"){
  idx <- which(adata$var$Symbol == gene)
  group <- adata$obs[, group]
  exprs_in <- t(adata$X)[idx,  ]
  dat <- data.frame("x" = as.factor(group), "y" = exprs_in)
  ggplot(dat, aes(x= x, y=y, fill = x)) + geom_violin(scale = "width") + 
    #ggbeeswarm::geom_quasirandom(size= 0.01)
    ylab("Expression") +
    xlab("") + theme_minimal() + theme(legend.position="none", axis.text.x = element_text(angle = 45, hjust = 1), axis.text = element_text(color = "black")) + ggtitle(gene)
}

#function to get average coordinates based on groups 
#' @param layout  - a layout to use, be it PCA, tSNE, UMAP, graph layout
#' @param groups - a factor to group by
average.coords <- function(layout, groups){
  groups <- as.factor(groups)
  dat <- data.frame("X" = layout[, 1], "Y" = layout[, 2])
  centroids <- do.call(rbind, lapply(levels(groups), function(x){
    cm <- colMeans(dat[groups == x, ])
    return(cm)
  }))
  rownames(centroids) <- levels(groups)
  return(centroids)
}

#annotated raster plot
#' @param layout the layout to use from the adata object
#' @param adata - the adata object to use
#' @param group - the grouping to use in adata$obs
annotated_raster_plot <- function(layout = "X_umap", group = "leiden", adata, xlab = "Dim1", ylab = "Dim2", size = 0.05, dpi = 600){
  require(ggrastr)
  group <- adata$obs[, group]
  layout <- adata$obsm$get(layout)
  av.coords  <- average.coords(layout, groups = group)
  set.seed(100)
  idx <- sample(1:nrow(layout))
  ggplot(data.frame("x" = layout[idx, 1], "y" = layout[idx, 2], "Factor" = as.factor(group[idx])), 
         aes(x = x, y=y, col = Factor)) + ggrastr::geom_point_rast(pch=19, size=size, raster.dpi = dpi) + xlab(xlab) +ylab(ylab)+ theme_classic() + 
    theme(legend.position="none") + annotate("text", x = av.coords[, 1], y=av.coords[, 2],label= rownames(av.coords))
}

#annotated raster plot
#' @param layout the layout to use from the adata object
#' @param adata - the adata object to use
#' @param group - the grouping to use in adata$obs
raster_plot <- function(layout = "X_umap", group = "leiden", adata, xlab = "Dim1", ylab = "Dim2", size = 0.05, dpi = 600){
  require(ggrastr)
  require(ggrastr)
  group <- adata$obs[, group]
  layout <- adata$obsm$get(layout)
  set.seed(100)
  idx <- sample(1:nrow(layout))
  ggplot(data.frame("x" = layout[idx, 1], "y" = layout[idx, 2], "Factor" = as.factor(group[idx])), 
         aes(x = x, y=y, col = Factor)) + ggrastr::geom_point_rast(pch=19, size=size, raster.dpi = dpi) + xlab(xlab) +ylab(ylab)+ theme_classic() + 
    guides(colour = guide_legend(override.aes = list(size=5)))
}

#plot HVG
#' @param adata - adaa object to use
hvg_plot <- function(adata){
  df <- adata$var
  ggplot(df, aes(x = means, y = dispersions_norm, color = highly_variable)) + ggrastr::geom_point_rast(size = 0.5) + 
    scale_color_manual(values = c("TRUE" = "black", "FALSE" = "grey50")) + theme_classic() + xlab("mean") + ylab("normalised dispersion") + 
    theme(legend.position = "none", axis.text = element_text(color = "black"))
}

#min max normalisation
#' @param x - a vector of values
min_max_normalisation <- function(x){
  return((x - min(x))/(max(x) -min(x)))
}

#a clean fractional heatmap function
#' @param adata - the adata object to use
fractional_heatmap <- function(adata, group = "leiden", genes, dendros = TRUE, min_max = FALSE, 
                               colorscheme = c("darkblue", "grey80", "red"),
                               filter_fraction = NULL, bias = 3){
  genes <-  unique(genes)
  exprs_in <- t(adata$X)
  rownames(exprs_in) <- adata$var$Symbol
  exprs_in <- exprs_in[!duplicated(rownames(exprs_in)), ]
  exprs_in <- exprs_in[rownames(exprs_in) %in% genes, ]
  
  #sort the clusters
  group <- adata$obs[, group]
  if(!is.factor(group)){group <- factor(group, levels = unique(group))}else{group <- factor(group, levels = levels(group))}
  #calculate fractions
  fractions <- do.call(cbind, lapply(levels(group), function(x){
    Matrix::rowSums(exprs_in[, group %in% x] > 0)/sum(group %in% x)
  }))
  colnames(fractions) <- levels(group)
  
  if(min_max){
    exprs_in <- t(apply(exprs_in, 1, min_max_normalisation))
    message("dotplot")
    message("min max normalised")
  }else{
    message("dotplot")
    message("not min max normalised")
  }
  #calculate averages
  averages <- do.call(cbind, lapply(levels(group), function(x){
    Matrix::rowMeans(exprs_in[, group %in% x])
  }))
  colnames(averages) <- levels(group)
  
    df <- stack(averages)
  df$fraction <-  stack(fractions)$value
  colnames(df) <- c("Gene", "Cluster", "Mean", "Fraction")
  if(dendros){ 
    rowclust <- hclust(dist(averages), method = "ward.D2")
    df$Gene <- factor(df$Gene, levels = rownames(averages)[rowclust$order])
    colclust <- hclust(dist(t(averages)), method = "ward.D2")
    df$Cluster <- factor(df$Cluster, levels = colnames(averages)[colclust$order])
    message("ordered by dendrograms")
  }else{
    df$Gene <- factor(df$Gene, levels = rev(genes)) 
    message("not ordered by dendrograms")
  }
  df$Fraction[df$Fraction < filter_fraction] <- NA
  
  #set the colorscheme
  colorscheme <- colorRampPalette(colors = colorscheme, bias = bias)(100)
  if(min_max){
    col_gradient <- scale_color_gradientn(colors = colorscheme, limits = c(0, 1))
  }else{
    col_gradient <- scale_color_gradientn(colors = colorscheme)
  }

  #plot it
 ggplot(data.frame(df), aes(x = Cluster, y = Gene, color = Mean, size = Fraction)) + 
    geom_point() + theme_bw() + col_gradient + scale_size_continuous(limits = c(0, 1)) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
          axis.text.y = element_text(color = "black", face = "italic" ), axis.line = element_blank(),
          panel.border = element_rect(colour = "black", size = 0.5)) + 
    ylab("") + xlab("") 

  
}

#a ggplot heatmap for a pre-ordered expression matrix
#' @adata - an anndata object
#' @scale_cap - maximal scaled value. will scale as a default
#' @cell_order - the ordering vector of cells
#' @dendros - whether to order the columns and rows by dendrograms - this can be expensive as we have to calculate a distance matrix
#' @colors - set of colors to use
#' @return  the heatmap
single_cell_heatmap <- function(adata, genes, cell_order, dendros = FALSE, scale_cap = 2.5, 
                           colors = list("low" =  "darkblue",
                                         "mid" = "grey90",
                                         "high" = "darkred")){
  
  genes <-  unique(genes)
  exprs_in <- t(adata$X)
  rownames(exprs_in) <- adata$var$Symbol
  exprs_in <- exprs_in[!duplicated(rownames(exprs_in)), ]
  expression_matrix <- as.matrix(exprs_in[rownames(exprs_in) %in% genes, ])
#scale the data
  heat_in <- t(scale(t(expression_matrix)))
  heat_in[heat_in > scale_cap] <- scale_cap
  heat_in[heat_in < -scale_cap] <- -scale_cap
  if(dendros){
    cell_order = NULL
    rowclust <- hclust(dist(heat_in), method = "ward.D2")
    colclust <- hclust(dist(t(heat_in)), method = "ward.D2")
    heat_in <- heat_in[rowclust$order, colclust$order]
    
  }else{
    colnames(heat_in) <- adata$obs_names$values
    heat_in <- heat_in[, cell_order]
  }
  #reshape it
  heat_in_mlt <- reshape2::melt(heat_in)
  levels(heat_in_mlt$Var2) <- colnames(heat_in)
  heat_in_mlt$Var1 <- factor(heat_in_mlt$Var1, levels=  rev(genes))
  ggplot(heat_in_mlt, aes(x = Var2, y = Var1, fill = value)) + ggrastr::geom_tile_rast() +
    scale_fill_gradient2(low = colors["low"],
                         mid = colors["mid"],
                         high = colors["high"],
                         midpoint = 0,
                         limit = c(-scale_cap,scale_cap)) + theme_classic() +
    theme(axis.text.y = element_text(color = "black", face = "italic", size = 7.5),  
          axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
          axis.ticks.y = element_line(colour = "black", size = 0.5),
          panel.border = element_blank(),
          axis.line = element_blank()
    ) + xlab("") + ylab("")
}

#tfidf markers

#' Calculates the tf-idf for a set of target cells against a background of the cells given in "universe".
#' From Matt Youngs script for kidney cells.
#'
#' @param data The data matrix to use.
#' @param target Columns that are the target.
#' @param universe Columns that we should consider (target must be a subset).
tfidf = function(data,target,universe){
  require(Matrix)
  if(!all(target %in% universe))
    stop('Target must be a subset of universe')
  nObs = Matrix::rowSums(data[,colnames(data) %in% target]>0)
  nTot = Matrix::rowSums(data[,colnames(data) %in% universe]>0)
  tf = nObs/length(target)
  idf = log(length(universe)/nTot)
  score = tf*idf
  #Calculate p-value for significance based on using a hypergeometric distribution to simulate the results of infinite random sampling
  pvals = phyper(nObs-1,nTot,length(universe)-nTot,length(target),lower.tail=FALSE)
  qvals = p.adjust(pvals,method='BH')
  ntf = (exp(-idf)*length(universe)-tf*length(target))/(length(universe)-length(target))
  return(data.frame(geneFrequency=tf,
                    geneFrequencyOutsideCluster=ntf,
                    geneFrequencyGlobal=exp(-idf),
                    geneExpression=Matrix::rowMeans(data[,colnames(data) %in% target]),
                    geneExpressionOutsideCluster = Matrix::rowMeans(data[, !colnames(data) %in% target]),
                    geneExpressionGlobal = Matrix::rowMeans(data),
                    idf=idf,
                    tfidf=score,
                    qval=qvals)[order(score,decreasing=TRUE),])
}



#find all markers with tfid
#' @param sce - an SCE object
#' @param groups - a set of cluster identities length of ncol(sce)
tfidf_all_markers <- function(adata, groups = "leiden", gene_subset = NULL, use_raw = FALSE){
  groups <- adata$obs[, groups]
  groups <- factor(groups)
  
  if(use_raw){
    adata <- adata$raw
  }else{adata <-adata}
  if(is.null(gene_subset)){
    gene_subset <- adata$var_names$values
    message("using all genes")
  }else{
    message("subsetting some genes")
  }
  #matrix organisation
  np <- import("numpy")
  mat_in <- np$transpose(adata$X)
  cellnames <- adata$obs_names$values
  colnames(mat_in) <-  cellnames
  rownames(mat_in) <- adata$var_names$values
  #then subset
  mat_in <- mat_in[rownames(mat_in) %in% gene_subset, ]
  library(pbmcapply)
  tfid.list <-  pbmclapply(levels(groups), function(i){
    tfid <- tfidf(data = mat_in, target = colnames(mat_in)[groups == i], universe = colnames(mat_in))
    tfid$gene <- adata$var[rownames(tfid), "Symbol"]
    return(tfid)
  }, mc.cores = 8)
  names(tfid.list) <- as.character(levels(groups))
  return(tfid.list)    
}


#to convert mouse gene lists
#' @param x - a vector of gene ids
#' @param type -  the attribute type
convertMouseGeneList <- function(x, type = "ensembl_gene_id"){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c(type), filters = type, values = x , mart = mouse, attributesL = c(type), martL = human, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  return(humanx)
}


mean_mean_plot <- function(adata, groups = "leiden", gp1 = "0", gp2 = "1", top = 10){
  #calculate means and difference in means
gp1_mns <- Matrix::colMeans(adata$X[adata$obs[, groups] %in% gp1,] )
gp2_mns <- Matrix::colMeans(adata$X[adata$obs[, groups] %in% gp2,])
df <- data.frame(gp1_mns, gp2_mns, diff = gp1_mns-gp2_mns, names = adata$var$Symbol)
df$label <- NA
df <- df[order(df$diff), ]
idx_labels <- c(1:top, (nrow(df)-top):nrow(df))
df$label[idx_labels] <-as.character( df$names[idx_labels])
#plot it
ggplot(df, aes(x = gp1_mns, y = gp2_mns, label = label, 
               fill = is.na(label))) + geom_point(pch = 21) + ggrepel::geom_label_repel(fill = "white", size =2 ) + 
  theme_bw() + scale_fill_manual(values = c("TRUE" = "grey50", "FALSE" = "red")) + 
  theme(legend.position = "none") + xlab(gp1) + ylab(gp2)
}

#this will plot acordinat to density
density_plot_input <- function(adata, use_rep = "X_pca", cell_variable){
  embedding <- adata$obsm$get(use_rep)
    df <- do.call(rbind, 
                  lapply(unique(adata$obs[, cell_variable]), function(x){
                    embed <- embedding[adata$obs[, cell_variable] %in% x, 1:2]
                    dens <- get_density(x = embed[, 1],y= embed[, 2])
                    return(data.frame(x= embed[, 1], y =embed[, 2], density = dens, variable = x ))
                  }))
  return(df)}
plot_group_density <- function(density_input){
  ggplot(density_input, aes(x=x,y=y, color = density)) + geom_point_rast(size = 0.1) + scale_color_viridis_c() + facet_wrap(.~variable) + 
    coord_fixed() + theme_void() + xlab("PC1") + ylab("PC2")
}

#simple clustering
#simple clustering beginning with a table of counts
#' @toc - the table of counts
#' @return - a vector of clusters
simple_pp_sc_clustering <- function(adata){
  cl_adata <- adata$copy()
  cl_adata <- cl_adata[cl_adata$obs$dd]
  sc$pp$normalize_per_cell(cl_adata)
  sc$pp$log1p(cl_adata)
  sc$pp$highly_variable_genes(cl_adata, n_top_genes = as.integer(3000))
  sc$pp$pca(cl_adata)
  sc$pp$neighbors(cl_adata, n_neighbors = as.integer(15))
  sc$tl$leiden(cl_adata, resolution = 0.5)
  annotated_raster_plot(adata= cl_adata, layout = "X_pca", group = "leiden")
  return(cl_adata$obs$leiden)
}

#do adata soupX
#a vanilla soupX pipeline to use when reading in the sample
#' @adata - the adata object, unfiltered but with default drops calculated
#' @return - a list of uncorrected and corrected adata object
do_soupx <- function(adata){
  np <- import("numpy")
  tod <- np$transpose(adata$X)
  cellnames <- adata$obs_names$values
  colnames(tod) <-  1:ncol(tod)
  rownames(tod) <- adata$var_names$values
  tod <- as(tod, "dgTMatrix")
  toc <- tod[, adata$obs$dd]
  soupchan = SoupChannel(tod,toc, calcSoupProfile = TRUE)
  #do a simple clustering
quick_clusters <- simple_pp_sc_clustering(adata)
  names(quick_clusters) <- colnames(toc)
  soupchan <- setClusters(soupchan, quick_clusters)
  soupchan = autoEstCont(soupchan, tfidfMin = 0.75, soupQuantile = 0.75, doPlot = FALSE)
  out = adjustCounts(soupchan, clusters = quick_clusters, verbose = 3,
                     roundToInt=TRUE)
  #subset to the cells that we are after
  adata <- adata[adata$obs$dd]
  #save the counts and corrected counts as seperate layers
  adata$layers <- list("counts" = adata$X, 
                       "corrected_counts" = t(out))
  return(adata)
}


#scrublet wrapper
#' @adata the anndata object
#' @channel the obs column that corresponds to the 10X channel. You could potentially use "donor" here...
scrublet_wrapper <- function(adata, channel = "channel"){
  scrb <- import("scrublet")
  scrb_results <- pblapply(unique(adata$obs[, channel]), function(x){
    adata_sub <-  adata[adata$obs[, channel] %in% x]
    scrub = scrb$Scrublet(adata_sub$X)
    scrub_out = scrub$scrub_doublets()
    return(scrub_out)
  })
  adata$obs$scrublet_scores <-  unlist(lapply(scrb_results, function(x){x[1]}))
  adata$obs$scrublet_cuts <-  unlist(lapply(scrb_results, function(x){x[2]}))
  return(adata)
}

#scrublet wrapper which takes into account a genotype doublet determined predicted doublet rate
#' @adata the anndata object
#' @channel the obs column that corresponds to the 10X channel. You could potentially use "donor" here...
scrublet_wrapper_w_genotypes <- function(adata, channel = "channel"){
  scrb <- import("scrublet")
  scrb_results <- pblapply(unique(adata$obs[, channel]), function(x){
    adata_sub <-  adata[adata$obs[, channel] %in% x]
    expected <- sum(!adata_sub$obs$status %in% "singlet")/length(adata_sub$obs$status) #program in an expected doublet rate which will be the proportion of cells which are genotype doublets
    if(expected == 0){
      expected = 0.05
    }else{expected = expected}
    scrub = scrb$Scrublet(adata_sub$X, expected_doublet_rate = expected)
    scrub_out = scrub$scrub_doublets()
    return(scrub_out)
  })
  adata$obs$scrublet_scores <-  unlist(lapply(scrb_results, function(x){x[1]}))
  adata$obs$scrublet_cuts <-  unlist(lapply(scrb_results, function(x){x[2]}))
  return(adata)
}

#test clusters for over-representation with doublets.
#scrublet_hypergeo
#' @scrublet_cuts - vector of TRUE FALSE boolean
#' @cluster - vector of cluster annotations
scrublet_test <- function(scrublet_cuts, clusters){
  tbl = table(clusters, scrublet_cuts) 
  fractional_tbl <- tbl/rowSums(tbl)
  q = tbl[, 2] #this is the number of TRUE in each cluster
  m = sum(scrublet_cuts) #this is the total number of TRUE in the dataset
  n = length(scrublet_cuts) - m #this is the total number of FALSE in the dataset
  k = table(clusters) # this is the number of cells in each cluster
  hypergeo_out <- phyper(q=q-1, m=m, n=n, k=k, lower.tail = FALSE)
  df <- data.frame("fraction_positive" = fractional_tbl[, 2],
             "pvalue" = hypergeo_out, 
             "p.adj"= p.adjust(hypergeo_out, method = "BH"),
             "cluster" = names(q), row.names =  names(q))
  return(df)
}


#do multibatch PCA
multibatch_pca <- function(adata, group){
  library(batchelor)
  library(BiocSingular)
  set.seed(100)
  mbpca_out <- batchelor::multiBatchPCA(t(adata$X)[adata$var$highly_variable, ], 
                                        batch = adata$obs[, group], 
                                        preserve.single = TRUE,
                                        BSPARAM = IrlbaParam(deferred = TRUE),
                                        BPPARAM = SerialParam())
  adata$obsm <- list("X_pca" = mbpca_out@listData[[1]])
  return(adata)
}

#for milo analysis
#fast countcells
fast_countCells <- function(x, samples, meta.data=NULL){
  
  # cast dplyr objects to data.frame
  if(!is.data.frame(meta.data) & !is.null(meta.data)){
    meta.data <- as.data.frame(meta.data)
  }
  
  if(length(samples) > 1 & !is.null(meta.data)){
    stop("Multiple sample columns provided, please specify a unique column name")
  } else if(is.null(meta.data) & length(samples) != ncol(x)){
    stop(paste0("Length of vector does not match dimensions of object. Length:",
                length(samples), " Dimensions: ", ncol(x)))}
  
  # check the nhoods slot is populated
  if(length(nhoods(x)) == 0){
    stop("No neighbourhoods found. Please run makeNhoods() first.")
  }
  
  message("Checking meta.data validity")
  if(!is.null(meta.data)){
    samp.ids <- unique(meta.data[, samples])
  } else{
    samp.ids <- unique(samples)
  }
  
  n.hoods <- length(nhoods(x))
  message(paste0("Setting up matrix with ", n.hoods, " neighbourhoods"))
  
  #precompute the cells in each sample
  sample_idx <- lapply(samp.ids, function(s){
    sample_idx <-  meta.data[, samples] %in% s
    return(c(1:nrow(meta.data))[sample_idx])
  })
  names(sample_idx) <- samp.ids
  
  message("Counting cells in neighbourhoods")
  nhood_by_sample <- pblapply(1:n.hoods, function(i){
    idx <- nhoods(x)[[i]]
    sample_counts <- list()
    for(j in samp.ids){
      sample_cells <- sample_idx[[j]]
      inters_l <- length(intersect(idx, sample_cells))
      sample_counts[[j]] <- inters_l
    }
    neigh_counts <- unlist(sample_counts)
    return(neigh_counts)
  })
  count.matrix <- do.call(rbind, nhood_by_sample)
  # add to the object
  rownames(count.matrix) <- c(1:n.hoods)
  nhoodCounts(x) <- count.matrix
  return(x)
} #this function improves the cell counting speed




#function to recluster existing leiden clusters on a graph.
recluster_leiden <- function(adata, clusters, resolution = 0.2){
  clusters <- as.list(clusters)
  sc$tl$leiden(adata, 
               restrict_to = 
                 tuple("leiden", clusters), 
               resolution = as.numeric(resolution))
  adata$obs$leiden <- as.character(adata$obs$leiden)
  adata$obs$leiden_R <- as.character(adata$obs$leiden_R)
  adata$obs$leiden[adata$obs$leiden %in% clusters ] <- adata$obs$leiden_R[adata$obs$leiden %in% clusters]
  adata$obs$leiden <- as.factor(adata$obs$leiden)
  return(adata)
}

