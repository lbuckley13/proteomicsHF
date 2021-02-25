linear.rf <- function(data, seed) {

  tune.tree <-
    randomForestSRC::tune(outc ~ ., as.data.frame(data),
                          mtryStart = 10, stepFactor = 1.25,
                          nodesizeTry = c(15:25, seq(30,300, by = 15)),
                          ntreeTry = 25, doBest = FALSE)

  rand.tree <-
    randomForestSRC::rfsrc(outc ~ ., as.data.frame(data),
                           nodesize = tune.tree$optimal[1],
                           mtry = tune.tree$optimal[2],
                           ntree = 500, block.size = 1,
                           seed = seed)

  list(rand.tree = rand.tree,
       tune.tree = tune.tree,
       vars.tree = colnames(data)) %>% return()
}


linear.filter.rf.array <- function(data, outcomes, iter, adjust, filt = 0.8) {

  list.of.tree.lists <-
    tibble::tibble(outcome = outcomes) %>%
    dplyr::mutate(forest.list = purrr::map(outcome, ~ linear.filter.rf(data, .x, iter, adjust, filt)),
                  forest.scre = purrr::map2(outcome, forest.list, ~ linear.rf.tree.score(data, .y, .x)))

  return(list.of.tree.lists)
}


linear.filter.rf <- function(data, outc, iter, adjust, filt = 0.8) {
  tictoc::tic()

  rand.seed <- 101010
  tree.list <- list()
  outcome <- outc
  rf.data <- data %>%
    dplyr::mutate(outc = data[[outc]]) %>%
    dplyr::select(-dplyr::all_of(outcome))

  for (i in 1:iter) {
    tree <- tree.list[[i]] <- rf(rf.data, rand.seed)
    vsel <- randomForestSRC::var.select(tree$rand.tree, verbose = FALSE)

    next.vars <- vsel$varselect %>%
      dplyr::mutate(vars = rownames(vsel$varselect)) %>%
      dplyr::filter(depth > quantile(depth, filt)) %>%
      dplyr::pull(vars)
    deletions <- next.vars[!(next.vars %in% adjust)]

    rf.data <- rf.data %>%
      dplyr::select(-dplyr::all_of(deletions))

  }

  tictoc::toc()
  return(tree.list)
}




linear.rf.tree.score <- function(data, min.tree, outc) {

  s.val <- min.tree$predicted.oob

  score <- log(s.val + 0.01) / sd(log(s.val + 0.01))
  outcome <- data[[outc]]
  lin.n <- lm(outcome ~ score)
  rsqred <- summary(lin.n)$adj.r.squared

  return(list(score = score, lin.n = lin.n, rsqred = rsqred))
}




linear.plot.var.imp <- function(min.tree, labels) {

  library(ggplot2)

  varselect <- randomForestSRC::var.select(min.tree, verbose = FALSE)
  plot_labs <- tibble::tibble(term = rownames(varselect$varselect),
                              name = labels[term,'name']) %>%
    dplyr::mutate(name = dplyr::if_else(is.na(name), term, name)) %>%
    dplyr::pull(name) %>%
    substr(1, 35)
  names(plot_labs) <- rownames(varselect$varselect)

  p <- plot(ggRandomForests::gg_minimal_depth(varselect), lbls = c(plot_labs)) +
    ggtitle("Minimal Depth Variable Selection Refined")

  return(p)
}
