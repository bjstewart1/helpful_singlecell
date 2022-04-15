get_cpdb_df = function(cpdb_dir, celltypes = NULL){
    cpdb_files <- lapply(list.files(cpdb_dir), function(f){
    read.table(file.path(cpdb_dir, f), sep = "\t", quote = NULL, header = TRUE)
  })  
    names(cpdb_files) <- gsub(".txt", "", list.files(cpdb_dir))
        #organise these into index, means, pval
    cpdb_index <- data.frame(cpdb_files$means[, 1:11],row.names = cpdb_files$means[, 1])
    means <- reshape::melt(data.frame(cpdb_files$means[, 12:ncol(cpdb_files$means)], row.names = cpdb_index$id_cp_interaction, "interaction" = cpdb_index$id_cp_interaction))
    p_vals <- reshape::melt(data.frame(cpdb_files$pvalues[, 12:ncol(cpdb_files$pvalues)], row.names = cpdb_index$id_cp_interaction, "interaction" = cpdb_index$id_cp_interaction))
    #make an overall interaction dataframe
    interaction_df = data.frame('interaction' = means$interaction, 'cells' = means$variable,
                            'mean' = means$value, 'pval' = p_vals$value)
    #subset to significant interactions
    
    #get the possible combinations of cell types
    combinations = list()
    p=1
    for(i in celltypes){
        for(j in celltypes){
            i = gsub("-", ".", i)
            j = gsub("-", ".", j)
        combinations[[p]] = data.frame("ct1" = i,  "ct2" = j)
                p = p+1
    }}
combinations = do.call(rbind, combinations)
combinations$combination = paste0(combinations$ct1, ".", combinations$ct2)
combinations$ct1 = gsub("[.]", "-", combinations$ct1 )
combinations$ct2 = gsub("[.]", "-", combinations$ct2 )
    
    #remove entirely insignificant interactions 
    message("removing entirely insignificant interactions")
    insignificant_interactions = unique(interaction_df$interaction)[unlist(lapply(unique(interaction_df$interaction), function(x){
    sum(interaction_df[interaction_df$interaction %in% x, 'pval']) == sum(interaction_df$interaction %in% x)}))]     
    interaction_df = interaction_df[!interaction_df$interaction %in% insignificant_interactions, ]

    message("adding cell-cell info")
    #write in the combinations to interaction df
    idx = match(interaction_df$cells, combinations$combination)
    interaction_df[, 'sender'] = combinations[idx, 1]
    interaction_df[, 'receiver'] = combinations[idx, 2]
    interaction_df$cells = paste0(interaction_df$sender, "|", interaction_df$receiver)    

        
    #add the interactions 
    message("adding interactions/molecules")
    cpdb_index <- data.frame(cpdb_files$means[, 1:11],row.names = cpdb_files$means[, 1])
    interactions <- cpdb_index[interaction_df$interaction, ]
    
    #then get the molecules
    interaction_df$molecule_a <- unlist(lapply(strsplit(as.character(interactions$interacting_pair), "[_]"), function(x){x[1]}))
    interaction_df$molecule_b <- unlist(lapply(strsplit(as.character(interactions$interacting_pair), "[_]"), function(x){x[2]}))

    cpdb_df = cbind(interaction_df, interactions)
    
    
    
#now we want to get the ligands and receptors the right way round because inexplicably cpdb get's this all askew
  cpdb_df$receptor_cell <- NA
  cpdb_df$receptor_cell[cpdb_df$receptor_a %in% "True"] <- cpdb_df$sender[cpdb_df$receptor_a %in% "True"]
  cpdb_df$receptor_cell[cpdb_df$receptor_b %in% "True"] <- cpdb_df$receiver[cpdb_df$receptor_b %in% "True"]
  cpdb_df$ligand_cell <- NA
  cpdb_df$ligand_cell[cpdb_df$receptor_a %in% "False"] <- cpdb_df$sender[cpdb_df$receptor_a %in% "False"]
  cpdb_df$ligand_cell[cpdb_df$receptor_b %in% "False"] <- cpdb_df$receiver[cpdb_df$receptor_b %in% "False"]
  cpdb_df$receptor <- NA
  cpdb_df$receptor[cpdb_df$receptor_a %in% "True"] <- cpdb_df$molecule_a[cpdb_df$receptor_a %in% "True"]
  cpdb_df$receptor[cpdb_df$receptor_b %in% "True"] <- cpdb_df$molecule_b[cpdb_df$receptor_b %in% "True"]
  cpdb_df$ligand <- NA
  cpdb_df$ligand[cpdb_df$receptor_a %in% "False"] <- cpdb_df$molecule_a[cpdb_df$receptor_a %in% "False"]
  cpdb_df$ligand[cpdb_df$receptor_b %in% "False"] <- cpdb_df$molecule_b[cpdb_df$receptor_b %in% "False"]
  
  cpdb_df[cpdb_df$receptor_a == cpdb_df$receptor_b, "ligand_cell"] <- cpdb_df[cpdb_df$receptor_a == cpdb_df$receptor_b, "sender"]
  cpdb_df[cpdb_df$receptor_a == cpdb_df$receptor_b, "receptor_cell"] <- cpdb_df[cpdb_df$receptor_a == cpdb_df$receptor_b, "receiver"]
  cpdb_df[cpdb_df$receptor_a == cpdb_df$receptor_b, "ligand"] <- cpdb_df[cpdb_df$receptor_a == cpdb_df$receptor_b, "molecule_a"]
  cpdb_df[cpdb_df$receptor_a == cpdb_df$receptor_b, "receptor"] <- cpdb_df[cpdb_df$receptor_a == cpdb_df$receptor_b, "molecule_b"]
  
  cpdb_df$LR <- paste0(cpdb_df$ligand, "->", cpdb_df$receptor)
  cpdb_df$SR <- paste0(cpdb_df$ligand_cell, "->", cpdb_df$receptor_cell)
    return(cpdb_df)
     }
 

#plot interaction counts...
plot_interaction_counts = function(interaction_df, celltypes, p_val_cutoff = 0.05){
interaction_df = interaction_df[interaction_df$pval < p_val_cutoff, ]

#make a dataframe of counts of interactions
counts_df = data.frame(table(interaction_df$cells))
#acquire senders and receivers
counts_df$sender = unlist(lapply(strsplit(as.character(counts_df$Var1),  "[|]"), function(x){x[1]}))
counts_df$receiver = unlist(lapply(strsplit(as.character(counts_df$Var1),  "[|]"), function(x){x[2]}))

#write these into a matrix of all possible combinations
matrix = matrix(0, ncol = length(celltypes), nrow = length(celltypes))
colnames(matrix) = rownames(matrix) = rev(celltypes)

#make a matrix of frequencies including zero values
for(i in rownames(counts_df)){
    sender = counts_df[i, "sender"]
    receiver = counts_df[i, "receiver"]
    matrix[rownames(matrix) %in% sender, colnames(matrix) %in% receiver] = as.numeric(counts_df[i, 'Freq'])
}
#then get this into long form
mat_plot = reshape2::melt(matrix)
    
#then get this into long form
mat_plot = reshape2::melt(matrix)
mat_plot$Var1 = factor(mat_plot$Var1, levels = celltypes)
mat_plot$Var2 = factor(mat_plot$Var2, levels =  celltypes)
  
#then plot it.
ggplot(mat_plot, aes(x = Var1, y = Var2, fill = value)) + geom_tile() + 
    scale_fill_gradientn(colors = viridis::magma(10)) + coord_fixed() + 
    theme(axis.text.x = element_text(color ='black', angle = 90, hjust = 1), 
      axis.text.y = element_text(color = 'black')) + xlab("sender") + ylab("receiver")

}


#plot interaction counts...
plot_interaction_counts_symmetric = function(interaction_df, celltypes, p_val_cutoff = 0.05){
#get an interaction df
interaction_df = interaction_df[interaction_df$pval < p_val_cutoff, ]

freq_table = data.frame(table(interaction_df$SR))
freq_table[, 'cell1'] = unlist(lapply(strsplit(as.character(freq_table$Var1), "->"), function(x){x[1]}))
freq_table[, 'cell2'] = unlist(lapply(strsplit(as.character(freq_table$Var1), "->"), function(x){x[2]}))
#do this with igraph - this will collect the edges into a symmetric matrix.
gr = igraph::graph_from_edgelist(as.matrix(freq_table[, c('cell1', 'cell2')]), directed=FALSE)
E(gr)$weight = freq_table$Freq
adjmat = igraph::as_adjacency_matrix(gr, attr = 'weight')   
if(isSymmetric(as.matrix(adjmat))){message("graph is symmetric")}
adjmat = as.matrix(adjmat) #convert from sparse
mat_plot = reshape2::melt(adjmat)
mat_plot$Var1 = factor(mat_plot$Var1, levels = celltypes)
mat_plot$Var2 = factor(mat_plot$Var2, levels =  celltypes)
ggplot(mat_plot, aes(x = Var1, y = Var2, fill = value)) + geom_tile() + 
    scale_fill_gradientn(colors = viridis::magma(10)) + coord_fixed() + 
    theme(axis.text.x = element_text(color ='black', angle = 90, hjust = 1), 
      axis.text.y = element_text(color = 'black'))
}


#cellphone DB plot
cellphone_db_plot <- function(cpdb_df, celltypes, senders, receivers, cluster_interactions = TRUE, 
                              filter_secreted  = FALSE, scaling = c("min_max", "scale"), remove_auto = TRUE,
                              filter_integrins = TRUE, 
                              annotations_drop = c("guidetopharmacology.org")){
    
      #get these ordered by how we provided them
  egrid <- expand.grid(receivers, senders)
  cpdb_df$SR <- factor(cpdb_df$SR, levels =  paste0(egrid$Var2, "->", egrid$Var1))
  
  #subset to just secreted interactions
  if(filter_secreted){
    cpdb_df <- cpdb_df[cpdb_df$secreted %in% "True", ]
  }
  
  #remoove annotation strategy
  for(i in annotations_drop){
    cpdb_df <- cpdb_df[!cpdb_df$annotation_strategy %in% i, ]
  }
  
  #remove auto interactions
  if(remove_auto){
    cpdb_df <- cpdb_df[!cpdb_df$sender == cpdb_df$receiver, ]
    cpdb_df <- cpdb_df[!cpdb_df$molecule_a == cpdb_df$molecule_b, ]
  }
  
  if(filter_integrins){
    cpdb_df <- cpdb_df[cpdb_df$is_integrin %in% "False", ]
  }
  
    
  #subset to the cells we are interested in 
  cell_types <- union(senders, receivers)
  cpdb_df <- cpdb_df[cpdb_df$sender %in% cell_types + cpdb_df$receiver %in% cell_types == 2, ]
  
  #subset to senders and receivers
  cpdb_df <- cpdb_df[cpdb_df$ligand_cell %in% senders, ]
  cpdb_df <- cpdb_df[cpdb_df$receptor_cell %in% receivers, ]
  
  #remove totally nonsignificant interactions
  #now find the interacting pairs which are all non significant 
  message("filtering non significant interactions")
  nonsig_interactions <- unique(cpdb_df$interaction)[unlist(lapply(unique(cpdb_df$interaction), function(x){
    return(sum( cpdb_df[cpdb_df$interaction %in% x, "pval"]) == sum(cpdb_df$interaction %in% x))
  }))] #this gives us the p values that sum to 1
  
  cpdb_df <- cpdb_df[!cpdb_df$interaction %in% nonsig_interactions, ]
  #and render nonsignificant p values null
  cpdb_df$pval[cpdb_df$pval == 1] =NA
  
  #make sender receiver a factor
  cpdb_df$SR <- factor(cpdb_df$SR, levels = unique(cpdb_df$SR))
  
  #scale the data
  mat <- matrix(0, ncol = length(unique(cpdb_df$LR)), nrow = length(unique(cpdb_df$SR)))
  colnames(mat) <- unique(cpdb_df$LR)
  rownames(mat) <- unique(cpdb_df$SR)
  message("clustering and/or scaling step...")

    for(i in 1:nrow(cpdb_df)){
    lr = cpdb_df[i, 'LR']
    sr = cpdb_df[i, 'SR']
    mat[sr, lr] = cpdb_df[i, "mean"]
    }
    
  if(scaling == "min_max"){
    mat <- apply(mat, 2, min_max_normalisation)
  }
  if(scaling == "scale"){
    mat <- scale(mat)
  }
  melted_mat <- reshape2::melt(mat)
  rownames(melted_mat) <- paste0(melted_mat$Var1, melted_mat$Var2)
  cpdb_df$scaled <- melted_mat[paste0(cpdb_df$SR, cpdb_df$LR), "value"]  
  
  
  #now cluster the interactions
  if(cluster_interactions){
    message("clustering interactions")
    cpdb_df$LR <- factor(cpdb_df$LR, levels = colnames(mat)[hclust(dist(t(mat)), method = "ward.D")$order])
  }else{
    cpdb_df$LR <- factor(cpdb_df$LR)
  }
    
    
    #plot
  message("plotting result")
  if(scaling == "scale"){
    scaled_colors <- c("palegreen4", "grey80", "plum3")
    max_scaled <- max(abs(cpdb_df$scaled))
    pl <- ggplot(cpdb_df, aes(x = SR, y = LR, size = -log10(pval +1e-3 ), fill = scaled )) + geom_point(pch = 21) + theme_bw() +
      scale_size_continuous(limits = c(0, 4)) + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
            axis.text.y = element_text(color = "black"),
            axis.title = element_blank()) + scale_fill_gradient2(low = scaled_colors[1], mid = scaled_colors[2], high = scaled_colors[3],
                                                                 limits = c(-max_scaled, max_scaled ))
    
  }
  if(scaling == "min_max"){
    pl <- ggplot(cpdb_df, aes(x = SR, y = LR, size = -log10(pval +1e-3 ), fill = scaled)) + geom_point(pch = 21) + theme_bw() +
      scale_size_continuous(limits = c(0, 4)) + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
            axis.text.y = element_text(color = "black"),
            axis.title = element_blank()) + scale_fill_gradientn(colors = viridis::magma(100), limits = c(0, 1)) 
  }
  if(scaling == FALSE){
    pl <- ggplot(cpdb_df, aes(x = SR, y = LR, size = -log10(pval +1e-3 ), fill = log1p(mean))) + geom_point(pch = 21) + theme_bw() +
      scale_size_continuous(limits = c(0, 4)) + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
            axis.text.y = element_text(color = "black"),
            axis.title = element_blank()) + scale_fill_gradientn(colors = viridis::magma(100), limits = c(0, max(log1p(cpdb_df$mean)))) 
    
    
  }
  message("done")
  return(pl)

}


plot_interaction_graph = function(interaction_df, celltypes, p_val_cutoff, n_interactions_cutoff = 10){
interaction_df = interaction_df[interaction_df$pval < p_val_cutoff, ]
freq_table = data.frame(table(interaction_df$SR))
freq_table[, 'cell1'] = unlist(lapply(strsplit(as.character(freq_table$Var1), "->"), function(x){x[1]}))
freq_table[, 'cell2'] = unlist(lapply(strsplit(as.character(freq_table$Var1), "->"), function(x){x[2]}))
freq_table = freq_table[order(freq_table$Freq, decreasing = FALSE), ]
#do this with igraph - this will collect the edges into a symmetric matrix.
gr = igraph::graph_from_edgelist(as.matrix(freq_table[, c('cell1', 'cell2')]), directed=FALSE)
E(gr)$weight = freq_table$Freq
    E(gr)$weight = as.integer(E(gr)$weight)
gr = delete_edges(gr, E(gr)[E(gr)$weight < n_interactions_cutoff])
#order the vertices
gr <- permute(gr, match(V(gr)$name, celltypes))
#plot the graph
library(ggraph)
ggraph(gr, layout = 'linear', circular = TRUE) + 
    geom_edge_arc(aes(width = weight, color = weight), alpha = 0.5) + geom_node_point() +geom_node_label(aes(label = name)) + coord_fixed() + 
scale_edge_color_gradientn(colors = viridis::magma(10))
}