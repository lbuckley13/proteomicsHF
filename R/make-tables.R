make.tables <- function(data, fupt, diag, mods) {

  py <- survival::pyears(survival::Surv(data[[fupt]], data[[diag]]) ~ 1,
                         scale = 36525)

  table.templ <- readr::read_csv(mods) %>%
    dplyr::mutate(funlist = purrr::map(`Summary Function`,
                                       ~ switch(.x,
                                                "mean" = mean.sd,
                                                "median" = median.IQR,
                                                "percent" = num.percent)),
                  calls = purrr::invoke_map(funlist, Expression),
                  resul = purrr::map_chr(calls, ~ eval(parse(text = .x))))

  data[[fupt]] <- data[[fupt]] / 365.25
  med.fl <- parse(text = median.IQR(fupt)) %>% eval()
  tab <- tibble::tibble(`Table Var` = c("Event Rate", "N", "n", "Follow Up Time (Years)"),
                        `resul` = c(round(sum(data[[diag]]) / py$pyears, 2),
                                    length(data[[diag]]),
                                    sum(data[[diag]]),
                                    med.fl)) %>%
    dplyr::bind_rows(table.templ %>% dplyr::select(`Table Var`, resul))

  names(tab) <- c("Table 1.", "Patient Characteristics")
  return(tab)
}


median.IQR <- function(str) {
  exp <- paste0("glue::glue('{med[2]} [{med[1]}, {med[3]}]', med = round(quantile(data[['",
                str, "']], c(0.25, 0.5, 0.75)), 2))") %>% as.character() %>%
    return()
}

mean.sd <- function(str) {
  exp <- paste0("glue::glue('{mean.v} ± {std}', mean.v = round(mean(data[['",
                str, "']], 2)), std = round(sd(data[['", str, "']]), 2))") %>% as.character() %>%
    return()
}

num.percent <- function(str) {
  exp <- paste0("glue::glue('{num[2]} ({round(per, 2)}%)', num = summary(data[['",str,"']]),
                per = 100*num[2]/(num[1] + num[2]))") %>%
    as.character() %>%
    return()
}

