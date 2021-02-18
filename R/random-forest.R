#' Random Survival Forest Wrapper Function
#'
#' This function wraps the [randomForestSRC::rfsrc()] function. It first runs an
#' initial tuning algorithm, identifying the most appropriate parameters for
#' random forest model fit. Then, the model runs the more detailed random forest
#' tree generation step. The tree is then returned. The time variable must be
#' stored in a variable named time, and outcome in variable named outc. These
#' variable translations are done automatically in the functions accessed by
#' users.
#'
#' @param data Tibble of data containing the variables to be regressed.
#' @param seed The random seed to standardize random forest creation
#'
#' @return \describe{
#' \item{rand.tree}{The random forest grown in this algorithm}
#' \item{tune.tree}{The tuning forest and parameters derived}
#' \item{vars.tree}{{The variables used to grow this forest}}}
#' @export
#' @examples
#'
#' \dontrun{
#'
#' tree <- rf(visit.data, seed=101010)
#'
#' TRUE USAGE:
#'
#' tree.list <- filter.rf(visit.data, time="fuptime", outc="hfdiag",
#'                        iter = 30, adjust = adjustment.variables)
#' }
rf <- function(data, seed) {

  tune.tree <-
    randomForestSRC::tune(Surv(time, outc) ~ ., as.data.frame(data),
                          mtryStart = 10, stepFactor = 1.25,
                          nodesizeTry = c(15:25, seq(30,300, by = 15)),
                          ntreeTry = 25, doBest = FALSE)

  rand.tree <-
    randomForestSRC::rfsrc(Surv(time, outc) ~ ., as.data.frame(data),
                           nodesize = tune.tree$optimal[1],
                           mtry = tune.tree$optimal[2],
                           ntree = 500, block.size = 1,
                           seed = seed)

  list(rand.tree = rand.tree,
       tune.tree = tune.tree,
       vars.tree = colnames(data)) %>% return()
}

#' Validation Tree Creation
#'
#' This function generates a random forest from a set of derived important
#' variables. For example, if a set of important random forest variables were
#' calculated from the ARIC Visit 5 Data, this function then uses the HUNT data
#' to create a tree from those variables. This tree is then returned, evaluating
#' how these variables perform in informativeness in the new dataset.
#'
#' @param data The tibble containing the Validation dataset
#' @param time The time to event variable for this analysis
#' @param outc The event indicator variable for this analysis
#' @param valid.terms The pre-calculated important terms needed to generate the
#'   tree.
#'
#' @importFrom magrittr %>%
#'
#' @return \describe{ \item{tree}{This is the Random Forest grown using the
#' supplied terms} \item{vsel}{This is the calcuation of relative variable
#' importance, showing how the provided terms replicate in the validation
#' dataset.}}
#' @export
#'
#' @examples
#'
#' \dontrun{
#'
#' hunt.data <- haven::read_dta('hunt_data.dta')
#' time.var <- 'fuptime'
#' outc.var <- 'hfdiag'
#'
#' valid.terms <- readRDS('ARIC_random_forest_derived_terms.RDS')
#'
#' validation.tree <- valid.rf(data = hunt.data, time = time.var,
#'                             outc = outc.var, valid.terms = valid.terms)
#'
#' }
#'
valid.rf <- function(data, time, outc, valid.terms) {
  tictoc::tic()

  rand.seed <- 101010
  outcomes <- c(time, outc)

  rf.data <- data %>%
    dplyr::mutate(time = data[[time]],
                  outc = data[[outc]]) %>%
    dplyr::select(dplyr::all_of(valid.terms),
                 -dplyr::all_of(outcomes))

  tree <- rf(rf.data, rand.seed)
  vsel <- randomForestSRC::var.select(tree$rand.tree, verbose = FALSE)

  next.vars <- vsel$varselect %>%
    dplyr::mutate(vars = rownames(vsel$varselect))

  tictoc::toc()

  return(list(tree = tree,
              vsel = vsel))

}

filter.rf <- function(data, time, outc, iter, adjust, filt = 0.8) {
  tictoc::tic()

  rand.seed <- 101010
  tree.list <- list()
  outcomes <- c(time, outc)

  rf.data <- data %>%
    dplyr::mutate(time = data[[time]],
                  outc = data[[outc]]) %>%
    dplyr::select(-dplyr::all_of(outcomes))

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

select.min.tree <- function(tree.list) {

  errors <-
    purrr::map_dbl(tree.list, ~ .x$rand.tree$err.rate[.x$rand.tree$ntree])
  min.error <- tree.list[[which(errors == min(errors))]] %>%
    return()

}

predict.new.tree <- function(data, tree, time, outc) {

  outcomes <- c(time, outc)
  rf.data <- data %>%
    dplyr::mutate(time = data[[time]],
                  outc = data[[outc]]) %>%
    dplyr::select(-dplyr::all_of(outcomes)) %>%
    as.data.frame()

  pred <- predict(tree, newdata = rf.data,
                  block.size = 1, outcome = 'test')
  return(pred)
}


rf.tree.score <- function(data, min.tree, time, outc) {

  s.val <- min.tree$predicted.oob

  score <- log(s.val + 0.01) / sd(log(s.val + 0.01))
  outcome <- survival::Surv(data[[time]], data[[outc]])
  cox.n <- survival::coxph(outcome ~ score)
  cstat <- summary(cox.n)$concordance[[1]]

  return(list(score = score, cox.n = cox.n, cstat = cstat))
}


plot.var.imp <- function(min.tree, labels) {

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
