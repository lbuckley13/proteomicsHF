


lin.modl <- function(.data, outcome, protein, .adjust = adjust, extend = 100) {

  formula <- glue::glue("{outcome} ~ {protein} + {paste0(.adjust, collapse = '+')}")
  lin.mod <- lm(as.formula(formula), data = .data) %>%
    broom::tidy(conf.int = TRUE)

  estm <- lin.mod[2, c('estimate', 'conf.low', 'conf.high')] %>% unlist()
  pval <- lin.mod[2, 'p.value'] %>% unlist()

  LR <- formatC(estm, format = "f", 3, 5)
  pv <- formatC(pval, format = "e", 2)

  st <- floor(-log10(pval) + log10(0.05)) %>%
    {max(c(., -1) + 1)} %>%
    {strrep("*", .)} %>%
    substring(1, extend)

  ret <- list(outcome = outcome,
              estm = estm[1],
              c.low = estm[2],
              c.high = estm[3],
              pval = pval,
              desc = glue::glue("{LR[1]} ({LR[2]}, {LR[3]}) p={pv}{st}")) %>%
    return()

}


lin.arry <- function(proteins, data, lin.outc, adjust, .src = 1) {
  results <- tibble::tibble(term = proteins) %>%
    dplyr::mutate(model = furrr::future_map(term, ~ lin.modl(data, lin.outc, .x, adjust), .progress = TRUE)) %>%
    tidyr::unnest_wider(model) %>%
    dplyr::mutate(src = as.factor(.src))
  return(results)
}



lin.arry.aggr <- function(proteins, data, all.outc, adjust, labels, .src = 1) {

  tictoc::tic()
  final.results <-
    purrr::map(all.outc, ~ lin.arry(proteins, data, .x, adjust, .src)) %>%
    purrr::reduce(dplyr::bind_rows)

  final.results.num <- final.results %>%
    dplyr::select(term, outcome, estm, pval) %>%
    tidyr::pivot_wider(names_from = outcome,
                       names_glue = "{outcome}_{.value}",
                       values_from = c(estm, pval)) %>%
    dplyr::mutate(Name = labels[term, "name"]) %>%
    dplyr::select(sort(dplyr::current_vars())) %>%
    dplyr::select(term, Name, dplyr::everything())

  final.results.char <-  final.results %>%
    dplyr::select(term, outcome, desc) %>%
    tidyr::pivot_wider(names_from = outcome,
                       names_glue = "{outcome} Linear Regression",
                       values_from = desc) %>%
    dplyr::mutate(Name = labels[term, "name"]) %>%
    dplyr::select(sort(dplyr::current_vars())) %>%
    dplyr::select(term, Name, dplyr::everything())

  final.results <- final.results %>%
    dplyr::mutate(Name = labels[term, "name"]) %>%
    dplyr::select(sort(dplyr::current_vars())) %>%
    dplyr::select(term, Name, dplyr::everything())

  tictoc::toc()
  return(list(all.res   = final.results,
              numeric   = final.results.num,
              character = final.results.char))

}


lin.reg.terms.to.keep <- function(lin.reg.data, pvalue.filt) {

  terms.sign <- lin.reg.data %>%
    dplyr::filter(pval < 0.05 / pvalue.filt) %>%
    dplyr::group_by(outcome) %>%
    dplyr::arrange(pval) %>%
    dplyr::summarise(term)

  terms <- terms.sign %>%
    dplyr::summarise(candidates = paste0('"', term, '"', collapse = ", ")) %>%
    dplyr::mutate(cand = purrr::map(candidates, ~ paste0("c(", .x, ")")),
                  candidates = purrr::map(cand, ~ parse(text = .x) %>% eval())) %>%
    dplyr::select(outcome, candidates)

  return(terms)
}
