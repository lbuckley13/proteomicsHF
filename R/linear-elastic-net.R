linear.elastic.net <- function(data, terms, adjust, outc, .alpha = 1) {

  tictoc::tic()
  data.full <- data[complete.cases(data), ]
  outcome <- data.full[,outc] %>% dplyr::pull()

  mat <- data.full %>%
    dplyr::select(dplyr::all_of(c(terms, adjust))) %>%
    data.matrix()

  p.f <- (!colnames(mat) %in% adjust) %>%
    as.numeric()

  foldId <- seq(10) %>%
    rep(length.out = nrow(mat)) %>%
    sample()

  enet <- glmnet::cv.glmnet(x = mat, y = outcome, penalty.factor = p.f,
                            alpha = .alpha, parallel = TRUE, nfolds = 10,
                            foldid = foldId, trace = 0)

  coefs <- enet %>%
    coef(s = 'lambda.1se') %>%
    as.matrix()

  terms <-  coefs %>%
    tibble::as_tibble(rownames = 'covar') %>%
    dplyr::filter(`1` != 0, grepl('Seq', covar)) %>%
    dplyr::pull(covar)

  tictoc::toc()

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
    purrr::map(enets,
               ~ dplyr::filter(.x, as.character(terms) != "NULL")) %>%
    purrr::map(~ dplyr::mutate(.x, score.mod =
                 purrr::map2(coefs, outcome, ~ linear.get.score(data, .x, adjust, .y)))) %>%
    purrr::map(~ tidyr::unnest_wider(.x, score.mod)) %>%
    return()
}


linear.elastic.net.array <- function(data, adjust, outcomes, .alpha = 1, fil = 4) {

  results <- outcomes %>%
    dplyr::mutate(elasnet = purrr::map2(outcome, candidates, ~ linear.elastic.net(data, .y, adjust, .x, .alpha))) %>%
    tidyr::unnest_wider(elasnet)

  top.cand <-
    outcomes$candidates %>%
    purrr::flatten_chr() %>%
    table() %>%
    sort(decreasing = TRUE)

  full.cands <-
    names(top.cand)[which(top.cand > fil)]

  res.full.arr <- outcomes %>%
    dplyr::mutate(enet.all.cand = purrr::map(outcome, ~ linear.elastic.net(data, full.cands, adjust, .x, .alpha))) %>%
    tidyr::unnest_wider(enet.all.cand)

  list(indiv.enets = results,
       all.c.enets = res.full.arr) %>%
    return()

}
