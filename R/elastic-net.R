#' Main Elastic Net Wrapper Function
#'
#' This function is a wrapper function to the [glmnet::cv.glmnet()] function to
#' conduct LASSO / ElasticNet regularized regression. The regression optimizes
#' the C statistic, forces the adjustment variables, and runs the algorithm
#' using parallel processing.
#'
#' @param data The tibble containing the merged Protein and Study Data.
#' @param terms The terms to consider in the LASSO Regression.
#' @param adjust The vector of Adjustment variables to force
#' @param time A character vector of the Survival time to event outcome
#'   variable.
#' @param outc A character vector of the Survival event indicator outcome
#'   variable.
#' @param .alpha The optional alpha parameter to change LASSO to Elastic Net.
#'
#' @return
#' \describe{
#' \item{model}{The returned cv.glmnet LASSO Elastic Net Model.}
#' \item{coefs}{The coefficients of the Non-zero covariates retained.}
#' \item{terms}{The Final retained Proteins from this model.}
#' }
#' @export
#' @importFrom magrittr %>%
#' @examples
#' \dontrun{
#'
#'  fifth.lasso <- elastic.net(fifth.visit,
#'                             bonferroni$kept,
#'                             adjust, time, outc)
#'
#' }
elastic.net <- function(data, terms, adjust, time, outc, .alpha = 1) {

  tictoc::tic()
  surv <- survival::Surv(data[[time]], data[[outc]])

  mat <- data %>%
    dplyr::select(dplyr::all_of(c(terms, adjust))) %>%
    data.matrix()

  p.f <- (!colnames(mat) %in% adjust) %>%
    as.numeric()

  foldId <- seq(10) %>%
    rep(length.out = nrow(mat)) %>%
    sample()

  enet <- glmnet::cv.glmnet(mat, surv, family = "cox",
                            type.measure = 'C', penalty.factor = p.f,
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


#' LASSO Score Constructor
#'
#' This function uses the LASSO Model retained coefficients to recreate the
#' linear predictor. This predictor is then adjusted and regressed against the
#' outcome. This predictor is only based on the protein values, thus adjustment
#' variables must be included.
#'
#' @param data The tibble containing the merged proteomic and study data.
#' @param coefs The LASSO coefficients output from [elastic.net()].
#' @param adjust The vector of adjustment variables
#' @param time The character vector of the Survival outcome time to event
#'   variable
#' @param outc The character vector of the Survival outcome event indicator
#'
#' @return \describe{ \item{cox}{A Cox Model of the LASSO Score}
#' \item{score}{The actual double vector constructed score} }
#' @export
#' @importFrom magrittr %>%
#'
#' @examples
#' \dontrun{
#'
#' fifth.score <- get.score(fifth.visit, lasso$coefs, adjust,
#'                          "fuptime", "hfdiag")
#'
#'
#' }
get.score <- function(data, coefs, adjust, time, outc) {

  subj.data.numeric <- data %>%
    dplyr::select(dplyr::all_of(rownames(coefs))) %>%
    dplyr::mutate_if(is.factor, ~ as.numeric(.x) - 1)

  a.vars <- which(rownames(coefs) %in% adjust)

  multiplier <- t(coefs) %>%
    rep(each = dim(data)[1]) %>%
    matrix(nrow = dim(data)[1])

  score <- rowSums(subj.data.numeric * multiplier[, -a.vars])

  modeldf <- dplyr::select(data, dplyr::all_of(adjust))
  outcome <- survival::Surv(data[[time]], data[[outc]])
  cox.mod <- survival::coxph(outcome ~ score + ., data = modeldf)

  return(list(cox = cox.mod,
              score = score))

}
