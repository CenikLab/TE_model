load("human_TE_rho.rda")
load("annotations_12_4_22_human.rda")
rownames(annotations_12_4_22_human) = gsub("-", ".", rownames(annotations_12_4_22_human))
GO_list = read.csv("GO_id.csv")
human_term = annotations_12_4_22_human[, GO_list$GO_id]

ranking <- function(human_TE_rho, anno_df) {
  
  network = abs(human_TE_rho)
  ord = order(rownames(network))
  network = network[ord, ord]
  #missing=setdiff(rownames(anno_df_tmp),rownames(network))
  #missing_rows_df <- data.frame(matrix(rep(0, length(missing) * ncol(anno_df_tmp)), 
                                       #nrow = length(missing), ncol = ncol(anno_df_tmp)))
  #rownames(missing_rows_df) <- missing
  #colnames(missing_rows_df) <- colnames(anno_df_tmp)
  #anno_df=rbind(anno_df_tmp,missing_rows_df)
  match.lab = match(rownames(anno_df), rownames(network))
  filt.lab = !is.na(match.lab)
  filt.net = match.lab[filt.lab]
  network = network[filt.net, filt.net]
  
  df_fin = data.frame(matrix(ncol = nrow(network), nrow = ncol(anno_df)))
  colnames(df_fin) = rownames(network)
  rownames(df_fin) = colnames(anno_df)
  i=1
  for (GO in colnames(anno_df)) {
    print(i)
    anno = anno_df[, GO, drop = FALSE]
    genes.labels = as.matrix(anno[rownames(network), , drop = FALSE])
      for (gene in rownames(genes.labels[genes.labels==1,,drop=FALSE])) {
        df_fin[GO, gene] = "in the term"
      } 

  l <- dim(genes.labels)[2]
  g <- dim(genes.labels)[1]
  ab <- which(genes.labels != 0, arr.ind = TRUE)
  
  sumin <- crossprod(network, genes.labels)
  sumall <- matrix(apply(network, 2, sum), ncol = dim(sumin)[2], nrow = dim(sumin)[1])
  predicts <- sumin/sumall
  nans <- which(genes.labels == 1, arr.ind = TRUE)
  predicts[nans] <- NA
      
  predicts <- apply(abs(predicts), 2, rank, na.last = "keep", ties.method = "average")
  predicts_fin <- predicts/(g-dim(ab)[1])
  df_fin[GO,]=t(predicts_fin)
  i=i+1
  }
  return(df_fin)
}


GO_pred <- ranking(human_TE_rho, human_term)
save(GO_pred, file="human_GO_pred_rank_new.rda")


