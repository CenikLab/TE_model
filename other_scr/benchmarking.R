###change this association file
load("annotations_12_4_22_human.rda")
load("human_TE_rho.rda")

network=human_TE_rho
rownames(annotations_12_4_22_human)=gsub("-",".",rownames(annotations_12_4_22_human))
anno=annotations_12_4_22_human
genes.labels=as.matrix(anno)
network=abs(network)
ord <- order(rownames(network))
network <- network[ord, ord]
ord <- order(rownames(genes.labels))
genes.labels <- as.matrix(genes.labels[ord, ])
match.lab <- match(rownames(genes.labels), rownames(network))
filt.lab <- !is.na(match.lab)
filt.net <- match.lab[filt.lab]
network <- network[filt.net, filt.net]
genes.labels <- as.matrix(genes.labels[filt.lab, ])

l <- dim(genes.labels)[2]
g <- dim(genes.labels)[1]
ab <- which(genes.labels != 0, arr.ind = TRUE)
n <- length(ab[, 1])
###mask one gene each time
df_fin=data.frame(row.names = colnames(anno))
colnames(df_fin)=c("AUROC_mean","AUROC_median","AURIOC_uq","n_genes")
count=1
for (j in 1:l) {
  test.genes.labels <- matrix(nrow=g)
  d <- which(ab[, 2] == j)  # Which indices the genes are in this particular GO group
  t <- length(d)  # Total number of genes in the GO group
  test.genes.labels=matrix(genes.labels[,j], nrow = g, ncol = t)
  ab_tmp=which(genes.labels[,j,drop=FALSE] != 0, arr.ind = TRUE)
  for (i in 1:t) {
    c <- 1 + (i - 1)  # GO group to look at (ie column)
    test.genes.labels[ab_tmp[,'row'], c][i] <- 0
    #tmp[ab_tmp[d], c][i] <- 0
  }
  ###known filter
  #filter <- matrix(genes.labels, nrow = g, ncol = n * l)
  #filter <- matrix(nrow=g)
  d <- which(ab[, 2] == j)  # Which indices the genes are in this particular GO group
  t <- length(d)  # Total number of genes in the GO group
  filter=matrix(genes.labels[,j], nrow = g, ncol = t)

  sumin <- crossprod(network, test.genes.labels)

  sumall <- matrix(apply(network, 2, sum), ncol = dim(sumin)[2], nrow = dim(sumin)[1])

  predicts <- sumin/sumall

  nans <- which(test.genes.labels == 1, arr.ind = TRUE)
  predicts[nans] <- NA

  predicts <- apply(abs(predicts), 2, rank, na.last = "keep", ties.method = "average")
  
  negatives <- which(filter == 0, arr.ind = TRUE)
  positives <- which(filter == 1, arr.ind = TRUE)
  predicts[negatives] <- 0

  # print('Calculate ROC - np')
  np <- colSums(filter) - colSums(test.genes.labels)  # Postives
  # print('Calculate ROC - nn')
  nn <- dim(test.genes.labels)[1] - colSums(filter)  # Negatives
  # print('Calculate ROC - p')
  p <- apply(predicts, 2, sum, na.rm = TRUE)
  # print('Calculate ROC - rocN')
  rocN <- (p/np - (np + 1)/2)/nn
  
  tmp_mean <- mean(rocN)
  tmp_median <- median(rocN)
  tmp_upq <- quantile(rocN)[4]
  df_fin[j,"AUROC_mean"]=tmp_mean
  df_fin[j,"AUROC_median"]=tmp_median
  df_fin[j,"AURIOC_uq"]=tmp_upq
  df_fin[j,"n_genes"]=t
  print(count)
  count=count+1
  write.csv(df_fin,"auroc_te_GO_12_22_human.csv")
}

