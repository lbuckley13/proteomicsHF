linear.elastic.net <- function(data, terms, adjust, outc, lambda = "lambda.min", penalty) {

  tictoc::tic()

  message('\r', crayon::red(outc), appendLF = FALSE)
  flush.console()

  data.full <- data[complete.cases(data), ]
  outcome <- data.full[,outc] %>% dplyr::pull()

  mat <- data.full %>%
    dplyr::select(dplyr::all_of(c(terms, adjust))) %>%
    data.matrix()

  p.f <- (!colnames(mat) %in% penalty) %>%
    as.numeric()

  foldId <- seq(10) %>%
    rep(length.out = nrow(mat)) %>%
    sample()

  enet <- glmnet::cv.glmnet(x = mat, y = outcome, penalty.factor = p.f,
                            alpha = 1, parallel = TRUE, nfolds = 10,
                            foldid = foldId, trace = 0)

  coefs <- enet %>%
    coef(s = lambda) %>%
    as.matrix()

  terms <-  coefs %>%
    tibble::as_tibble(rownames = 'covar') %>%
    dplyr::filter(`1` != 0, grepl('Seq', covar)) %>%
    dplyr::pull(covar)

  tictoc::toc(quiet = TRUE)

  return(list(model = enet,
              coefs = coefs,
              terms = terms))
}



linear.get.score <- function(data, coefs, adjust, outc) {

  subj.data.numeric <- data %>%
    dplyr::select(dplyr::all_of(rownames(coefs)[-1])) %>%
    dplyr::mutate_if(is.factor, ~ as.numeric(.x) - 1)

  a.vars <- which(rownames(coefs) %in% adjust)

  multiplier <- t(coefs) %>%
    rep(each = dim(data)[1]) %>%
    matrix(nrow = dim(data)[1])

  score <- rowSums(subj.data.numeric * multiplier[, -a.vars])

  modeldf <- dplyr::select(data, dplyr::all_of(adjust))
  outcome <- data[[outc]]
  lin.mod <- lm(outcome ~ score + ., data = modeldf)
  rsq <- summary(lin.mod)$adj.r.squared

  return(list(linear = lin.mod,
              score = score,
              rsqred = rsq))
}


linear.get.score.enet <- function(enets, data, adjust) {

  results <-
    dplyr::mutate(enets, score.mod = purrr::map2(coefs, outcome, ~ linear.get.score(data, .x, adjust, .y))) %>%
    tidyr::unnest_wider(score.mod) %>%
    return()
}


linear.elastic.net.array <- function(data, adjust, outcomes,
                                     lambda = "lambda.min", force) {
  tictoc::tic()
  results <- outcomes %>%
    dplyr::mutate(elasnet = purrr::map2(outcome, candidates, ~ linear.elastic.net(data, .y, adjust, .x,
                                                                                  lambda = lambda,
                                                                                  penalty = force))) %>%
    tidyr::unnest_wider(elasnet)
  tictoc::toc()
  return(results)

}
