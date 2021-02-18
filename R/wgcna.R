soft.threshold <- function(data) {

  soft.threshold <-
    WGCNA::pickSoftThreshold(wg.data, powerVector = 1:20,
                             networkType = 'unsigned',
                             RsquaredCut = 0.9, verbose = 0)

  p <- threshold.plot(soft.threshold)

  return(list(sft = soft.threshold, p = p))
}

threshold.plot <- function(soft.threshold) {
  b <- soft.threshold$fitIndices[,1]

  p1 <- qplot(Power, SFT.R.sq, data = soft.threshold$fitIndices,
       xlab = "Soft Threshold (power)",
       ylab = "Scale Free Topology Model Fit, Unsigned R^2",
       main = "Scale Independence Model", geom="text", label = b) +
    geom_hline(yintercept = 0.9)

  p2 <- qplot(Power, truncated.R.sq, data = soft.threshold$fitIndices,
              xlab = "Soft Threshold (power)",
              ylab = "Scale Free Topology Model Fit, Unsigned R^2",
              main = "SSI Model More Complexity",
              geom="text", label = b) +
    geom_hline(yintercept = 0.9)

  p3 <- qplot(Power, mean.k., data = soft.threshold$fitIndices,
              xlab = "Soft Threshold (power)",
              ylab = "Mean Connectivity",
              main = "Mean Connectivity",
              geom="text", label = b)
  p <- gridExtra::grid.arrange(p1, p2, p3, ncol = 3)
  return(p)
}

hierTree <- function(data, sft) {

  tictoc::tic()

  pow <- sft$sft$powerEstimate

  adj <- WGCNA::adjacency(data, power = pow)
  dissTOM <- 1 - WGCNA::TOMsimilarity(adj)
  tree <- fastcluster::hclust(as.dist(dissTOM), method = "average")

  dynamicMods <- dynamicTreeCut::cutreeDynamic(dendro = tree,
                                               distM = dissTOM,
                                               deepSplit = 4,
                                               cutHeight = 0.9995,
                                               minClusterSize = 10,
                                               pamRespectsDendro = F)

  dynamicColors <- WGCNA::labels2colors(dynamicMods)

  tictoc::toc()

  WGCNA::plotDendroAndColors(tree, dynamicColors, "Dynamic Tree Cut",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      main = "Protein dendrogram and module colors")

  return(list(htree = tree,
              dynMd = dynamicMods,
              dynCl = dynamicColors))
}

summarize.modules <- function(data, colors, sft, labels) {

  tictoc::tic()

  pow <- sft$sft$powerEstimate

  hub.proteins <-
    WGCNA::chooseTopHubInEachModule(data, colors, omitColors = NA,
                                    power = pow, type = "unsigned")

  tictoc::toc()

  module.summary <- tibble::tibble(color = colors %>% table() %>% names(),
                                   freq  = colors %>% table(),
                                   hubs  = hub.proteins) %>%
    dplyr::mutate(Name = labels[hubs, 'name']) %>%
    return()
}

calculate.eigengenes <- function(data, colors) {

  tictoc::tic()

  MEList <- WGCNA::moduleEigengenes(data, colors = colors)
  eigenGenes <- tibble::tibble(MEList$eigengenes)

  MEDiss <- 1 - cor(eigenGenes)
  METree <- fastcluster::hclust(as.dist(MEDiss), method = "average")

  MEDissThres = 0.1

  plot(METree, main = "Module EigenGenes Clustering", xlab = "", sub = "")

  abline(h = MEDissThres, col = "red")

  merge <- WGCNA::mergeCloseModules(data, colors,
                             cutHeight = MEDissThres,
                             verbose = 3)
  tictoc::toc()

  return(list(eigenGenes = eigenGenes,
              mergedColors = merge$colors,
              mergedMEs = merge$newMEs))
}

plot.new.colors <- function(tree, colors1, colors2, labs) {

  WGCNA::plotDendroAndColors(tree, cbind(colors1, colors2),
                             labs, dendroLabels = FALSE, hang = 0.03,
                             addGuide = TRUE, guideHang = 0.05)

}
