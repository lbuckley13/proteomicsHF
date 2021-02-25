#' Cox Model Wrapper Function
#'
#' This function serves to wrap the survival coxph function for easier usage in
#' evaluating many cox models simultaneously. This function is accessed by
#' [cox.arry()] in order to parallelize the univariate comparisons of every
#' protein against the outcome.
#'
#' This function returns the effect estimate for the "protein" of interest
#' adjusted. The function also creates a string descriptor to print output as
#' table for \code{.csv} files. The significance of the p value is indicated by
#' each marker denoting a power of ten smaller.
#'
#' @param .data The data tibble containing the variables to be regressed.
#' @param time A character vector of the outcome variable time to event
#' @param outcome A character vector of the outcome variable event indicator
#' @param protein A character vector indicating which protein is being
#'   regressed.
#' @param .adjust A vector of all adjustment variables for this model.
#' @param extend Number of asterisks to add to Pvalue table indicator.
#'
#' @return A list with three components: \describe{ \item{hazr}{The Hazard Ratio
#'   of the Protein} \item{cint}{A Vector of length 2 with Lower and Upper
#'   Confidence Interval} \item{pval}{The P value of the model covariate}
#'   \item{desc}{A string that combines the relevant information (such as hazard
#'   ratio and confidence interval) needed to display in a table.}}
#' @export
#' @importFrom magrittr %>%
#'
#' @examples
#' \dontrun{
#'
#'   data <- haven::read_dta('~/proteomics/studydata.dta')
#'
#'   cox <- cox.modl(data, time = 'adjudhfdate', outcome = 'adjudhf_bwh',
#'                   protein = 'SeqId_7655_11', .adjust = adjust)
#'
#' }
#'
#'
cox.modl <- function(.data, time = "fuptime", outcome = "hfdiag",
                     protein, .adjust = adjust, extend = 100) {

  formula <- glue::glue("survival::Surv({time}, {outcome}) ~ {protein} + {paste0(.adjust, collapse = '+')}") %>%

  cox.mod <- survival::coxph(as.formula(fomula), data = .data) %>%
    broom::tidy(exponentiate = TRUE)

  hazr <- cox.mod[1, c('estimate', 'conf.low', 'conf.high')] %>% unlist()
  HR  <- formatC(hazr, format = "f", 3, 5)

  pval <- cox.mod[1, 'p.value'] %>% unlist()
  pv  <- formatC(pval, format = "e", 2)

  st <- floor(-log10(pval) + log10(0.05)) %>%
    {max(c(., -1) + 1)} %>%
    {strrep("*", .)} %>%
    substring(1, extend)

  ret <- list(hazr = hazr[1],
              cint = hazr[-1],
              pval = pval,
              desc = glue::glue("{HR[1]} ({HR[2]}, {HR[3]}) p={pv}{st}")) %>%
    return()
}


#' Cox Array Main Function
#'
#' This is the main access point for using parallel processing to complete
#' numerous Cox regressions simultaneously.
#'
#' This function utilizes the [furrr::future_map()] function to complete many
#' regressions quickly. Each protein is regressed individually against the time
#' and outcome variable specified, while adjusted for the provided variables.
#'
#' The .src input is specified to distinguish data between derivation data and
#' validation data. Derivation data should be \code{.src = 1} while validation
#' data should be \code{.src = 2}.
#'
#' @param proteins This is the vector of all proteins to be considered for
#'   Uniprotein, Adjusted Cox Regressions.
#' @param data This is the tibble dataframe with the covariates for regression
#' @param time This is a character vector of the Survival time to event Outcome.
#' @param outc This is a character vector of the Survival event indicator
#'   Outcome.
#' @param adjust This is a vector of all the adjustment variables.
#' @param .src This is an indicator of whether the models are the derivation or
#'   validation models. It is 1 when the data is the derivation dataset, and 2
#'   when it is the validation dataset.
#'
#' @return a Tibble containing the results of the multiple regression models.
#' @export
#' @importFrom magrittr %>%
#'
#' @examples
#' \dontrun{
#'   data <- haven::read_dta('~/proteomics/studydata.dta') %>% ...
#'       (merge study data with SomaData, see [get.visit()] for specifics)
#'
#'   univariate.results <- cox.arry(labels$term, data,
#'                                  time = "fuptime",
#'                                  outc = "hfdiag",
#'                                  adjust = adjust,
#'                                  .src = 1)
#'
#'  validation.results <- cox.arry(labels$term, data.validation,
#'                                  time = "fuptime",
#'                                  outc = "hfdiag",
#'                                  adjust = adjust,
#'                                  .src = 2)
#'
#' }
cox.arry <- function(proteins, data, time, outc, adjust, .src = 1) {
  tictoc::tic()
  results <- tibble::tibble(term = proteins) %>%
    dplyr::mutate(model = furrr::future_map(term, ~ cox.modl(data, time, outc, .x, adjust), .progress=TRUE)) %>%
    tidyr::unnest_wider(model) %>%
    dplyr::mutate(src = as.factor(.src))
  tictoc::toc()
  return(results)
}


#' Filtering Function to Identify which Proteins are retained
#'
#' This function calculates which proteins are retained at specified
#' significance intervals. For example, when the univariate results are
#' provided, this function returns all proteins that are bonferroni significant.
#'
#' This function can also be used to check which proteins are retained in
#' validation.
#'
#' @param data The Univariate (and Validation) Tibble with Regression Results.
#'   Output of [cox.arry()].
#' @param .pval The significance level to examine which proteins are retained.
#' @param .src The tibble in which we examine retained proteins. .src = 1 means
#'   the Derivation dataset if the derivation and validation sets are row-bound.
#' @param .terms Vector of terms to limit examination of significance to. Not
#'   necessary to be specified.
#'
#' @return
#' \describe{
#' \item{kept}{A vector containing all Proteins that were retained by the filter}
#' \item{data}{A Tibble with the Hazard Ratio and Significance of Retained Proteins}
#' }
#' @export
#' @importFrom magrittr %>%
#'
#' @examples
#' \dontrun{
#' plot.data <-
#'   dplyr::bind_rows(univariate.results,
#'                    validation.results) %>%
#'   dplyr::mutate(Name = labels[term, "name"])
#'
#' bonferroni <-
#'   terms.to.keep(plot.data,
#'                 .pval = 0.05/length(proteins), 1)
#' falsediscr <-
#'   terms.to.keep(plot.data,
#'                 .pval = 0.05, 1)
#'
#' retained <-
#'   terms.to.keep(plot.data,
#'                 .pval = 0.05, 2,
#'                 .terms = bonferroni$kept)
#'
#' }
terms.to.keep <- function(data, .pval, .src, .terms = NULL) {

  if (!is.null(.terms)) {
    data <- data %>% dplyr::filter(term %in% .terms)
  }

  kept <- data %>%
    dplyr::filter(pval <= .pval & src == .src) %>%
    dplyr::pull(term)

  dfrt <- data %>% dplyr::filter(term %in% kept)

  return(list(kept = kept, data = dfrt))
}


#' Volcano Plot Generator To Summarize Univariate Protein Significance
#'
#' This function takes in plot data and uses it to construct a Volcano Plot to
#' visualize significance. The legend is dependent on the derivation and
#' validation sets. As currently implemented, this function must have both
#' datasets bound.
#'
#' @param data The plot data (with vars \code{hazr} and \code{pval}) to create
#'   volcano plot with.
#' @param .pval1 The significance indicator for scaling 0.05 for the derivation
#'   dataset. The number of comparisons must be provided (i.e. 4877 for
#'   bonferroni or 1 for FDR).
#' @param .pval2 The significance indicator for scaling 0.05 for the validation
#'   dataset. The number of comparisons must be provided (i.e. 4877 for
#'   bonferroni or 1 for FDR).
#' @param suffix The Suffix that is used on the data join to separate the
#'   Derivation and Validation Sets.
#' @param ... Additional Parameters (such as title) to be passed to
#'   [EnhancedVolcano::EnhancedVolcano()]
#'
#' @return A Volcano Plot, comparing Hazard Ratio against -log10(pvalue)
#' @export
#' @importFrom magrittr %>%
#'
#' @examples
#' \dontrun{
#'
#' plot.data.visit.five <-
#'   dplyr::left_join(univariate.results,
#'                    validation.results,
#'                    by = "term",
#'                    suffix = c(".V5", ".V3"))
#'
#'   volcanoPlot(plot.data.visit.five,
#'               length(proteins),
#'               length(bonferroni$kept),
#'               suffix = c(".V5", ".V3"))
#'
#' }
volcanoPlot <- function(data, .pval1, .pval2, suffix, ...) {

  plot.df <- data
  prim.pv <- paste0("pval", suffix[1])
  secn.pv <- paste0("pval", suffix[2])

  keyvals <- plot.df[[prim.pv]] %>%
    cut(breaks = c(0, 0.05/.pval1, 0.05, 1),
        labels = c("royalblue", "red2", "grey30")) %>%
    as.character()

  names(keyvals) <- plot.df[[secn.pv]] %>%
    cut(breaks = c(0, 0.05/.pval1, 0.05, 1),
        labels = c("bfr-sig", "fdr-sig", "not-sig") %>%
          paste0(suffix[1]) %>%
          tolower()) %>%
    as.character()

  f <- plot.df[[secn.pv]] %>%
    cut(breaks = c(0, 0.05/.pval2, 0.05, 1),
        labels = c(17, 20, 4))

  shapevals <- as.numeric(levels(f))[f]

  names(shapevals) <- plot.df[[secn.pv]] %>%
    cut(breaks = c(0, 0.05/.pval2, 0.05, 1),
        labels = c("bfr-sig", "fdr-sig", "not-sig") %>%
          paste0(suffix[2]) %>%
          tolower()) %>%
    as.character()

  EnhancedVolcano::EnhancedVolcano(data, lab = data$term,
                                   x = paste0("hazr", suffix[1]), y = prim.pv,
                                   xlab = "Hazard Ratio", xlim = c(-0.2, 2.2),
                                   ylab = "-log10(p)",
                                   pCutoff = 0.05 / .pval1,
                                   FCcutoff = 0, pointSize = 3,
                                   selectLab = c("N"),
                                   hline = c(0.05, 0.05/.pval1, 0.05/.pval2),
                                   colCustom = keyvals,
                                   shapeCustom = shapevals,
                                   ...)
}
