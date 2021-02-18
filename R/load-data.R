#' Get Protein Names and Labels from SeqId Identifiers
#'
#' This is a function to generate a dictionary to translate from SeqId to
#' Protein Name and Entrez Gene Symbol. Output from this function is used in
#' creating tables, figures, and labels for SeqIds. This function also requires
#' an input for data cleaning flags. Proteins that are flagged will not be
#' included in the final label frame.
#'
#' @details The term, name, gene, and flag functions must be provided as-is (not
#' as strings). This is a quirk of using tidyverse and dplyr.
#'
#' @param path File Path to SomaLogic Data Dictionary
#' @param term Column Variable ID for the SeqId term
#' @param name Column Variable ID for the UniProt Full Name Column
#' @param gene Column Variable ID for the Entrez Gene Symbol Column
#' @param flag Filtering Function for Proteins that Pass QC
#'
#' @return Tibble with columns \code{term} for SeqId, \code{name} for
#'   UniProt Name, and \code{gene} for Entrez Gene Symbol. Data Frame rownames
#'   are set as the SeqId identifiers.
#'
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' \dontrun{
#' labels <- get.labels('data/proteinMapping.csv', sid,
#'                      uniprot.full.name,
#'                      entrezgenesymbol,
#'                      flag2 == 0)
#' }
#'
get.labels <- function(path, term, name, gene, flag) {
  data <-
    readr::read_csv(path) %>%
    dplyr::mutate(term = as.character({{ term }}),
                  temp = as.character({{ term }}),
                  name = as.character({{ name }}),
                  gene = as.character({{ gene }})) %>%
    dplyr::filter({{ flag }})

  tibble::column_to_rownames(data, var = "temp") %>%    # Translates SeqId into
    dplyr::select(term, name, gene) %>%                 # rownames for dictionary
    return()                                            # functionality

}

#' Load Raw SomaLogic Proteomic Data From Text File
#'
#' This is the function to read and process raw text file proteomics data from
#' SomaLogic. Data can provided with a label to preserve mapping from SeqId to
#' UniProt Full Name. Study ID is preserved in the ID data column.
#'
#' @param path File path to the SomaLogic text file with protein data.
#' @param labels Output of [get.labels()] function protein mapping dictionary
#' @param id Variable name of Study Identifier in SomaLogic Text File used to
#'   map study subjects to clinical data.
#'
#' @return Tibble containing raw proteomic data for all measured subjects.
#'   Tibble will be labelled if labelling dictionary is provided. Sample ID to
#'   map to subject study id is stored in the \code{id} variable.
#' @export
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#' labels <- get.labels('data/proteinMapping.csv', sid,
#'                      uniprot.full.name,
#'                      entrezgenesymbol,
#'                      flag2 == 0)
#'
#' soma.data.fifth <- get.soma('data/somaV5.txt', labels = labels)
#' soma.data.third <- get.soma('data/somaV3.txt', labels = labels)
#' }
get.soma <- function(path, labels = NULL, id = "SampleId") {
  data <-
    readr::read_tsv(path) %>%
    dplyr::mutate(id = .data[[id]]) %>%
    dplyr::select(dplyr::all_of(c('id', labels$term)))

  if(!is.null(labels)) {
    labelled::var_label(data) <- labels[names(data), 'name']
  }

  return(data)
}

#' Filter Study Data Set to Variables of Interest
#'
#' This function filters and creates new variables from the study data. This
#' function takes csv file inputs that detail the new variable names and
#' formulas created from the variable names in the Stata file.
#'
#' @details This function uses csv files to inject commands into
#'   [dplyr::mutate()] and [dplyr::filter()] functions. The format of the csv
#'   file to properly execute the desired commands is detailed below. The
#'   \code{mods} parameter is required and must be a list of path(s) to the csv
#'   files. Convention is to include one file with the outcome
#'   variables--including heart failure diagnosis variables and time to event
#'   variables and another file including adjustment variables, such as
#'   demographics and clinical history. The filter csv should include filtering
#'   conditions, such as excluding patients with prevalent heart failure, and
#'   the format of that csv is also detailed below.
#'
#' @section CSV Formatting:
#' ## Adjustment and Outcome Variables
#' Files must include a header for three columns--Variable, Expression, Label.
#' Then, a new row must be created for each new variable, its expression, and
#' its label. Expressions are R code snippets that calculate variables. Data in
#' tibble should be assumed to be attached--same as when tidyverse formats its
#' variables.
#'
#' !! First row in adjustment variables must always be the study
#' identifier for merging purposes. !!
#'
#' ## Examples:
#' ### Adjustment Variables:
#'
#' | Variable | Expression           | Label                     |
#' |----------|----------------------|---------------------------|
#' | id       | id                   | ARIC COHORT STUDY ID      |
#' | age      | v5age51              | Visit 5 Age               |
#' | bmi      | v5_bmi51             | Visit 5 BMI               |
#' | race     | as.factor(race == 1) | Subject Race (Black == 1) |
#'
#' This table is coded in csv as:\cr
#'
#' Variable, Expression, Label\cr
#' id, id, ARIC COHORT STUDY ID\cr
#' age, v5age51, Visit 5 Age\cr
#' bmi, v5_bmi51, Visit 5 BMI,\cr
#' race, as.factor(race == 1), Subject Race (Black == 1)\cr
#'
#' ### Outcomes Variables:
#'
#' | Variable | Expression                        | Label                     |
#' |----------|-----------------------------------|---------------------------|
#' | hfdiag   | !is.na(adjudhf_bwh)               | Incident Heart Failure Dx |
#' | fuptime  | as.double(adjudhfdate - v5date51) | V5 HF Follow Up Time      |
#'
#' This table is coded in csv as:\cr
#'
#' Variable, Expression, Label\cr
#' hfdiag, !is.na(adjudhf_bwh), Incident Heart Failure Dx\cr
#' fuptime, as.double(adjudhfdate - v5date51), V5 HF Follow Up Time\cr
#'
#' !! The First Row Must be the primary outcome, and the last row must be the
#' primary time to event variable. This quirk may be patched in the coming
#' versions.
#'
#' ## Filters CSV Format:
#' These files filter data according to exclusion criteria. Most used to remove
#' subjects with prevalent heart failure.
#'
#' ## Data Filters and Exclusions:
#' Files must include headers for three columns--Variable, Operator, Expression.
#' These headers represent the new variable, its filtering operation
#' (==, >, >=, <, <=, etc.), and the expression for this operation respectively.
#'
#' ## Examples:
#'
#' | Variable    | Operator | Expression |
#' |-------------|----------|------------|
#' | v5_prevhf52 | ==       | FALSE      |
#' | fuptime     | >        | 0          |
#'
#' This table is coded in csv as:\cr
#'
#' Variable, Operator, Expression\cr
#' v5_prevhf52, ==, FALSE\cr
#' fuptime, >, 0\cr
#'
#' @param stata Study (ARIC or HUNT) dataset in stata format
#' @param mods List of paths to formatted csv files creating new adjustment
#'   variables. Format of csv is detailed below. Files should be titled
#'   \code{adjusted.csv} and the other \code{outcomes.csv}. This parameter is
#'   required
#' @param filts List of path (singular) to formatted csv of the data filters to
#'   be applied to the dataset
#'
#' @return Tibble with data from all variables found in the mods csv files.
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' \dontrun{
#' fifth.visit.study <-
#'   get.adjust(stata = haven::read_dta('data/ARICmaster_121820.dta'),
#'              mods  = list('data/visit-five/outcomes.csv',
#'                           'data/visit-five/adjusted.csv'),
#'              filts = list('data/visit-five/filters.csv'))
#'
#' fifth.visit.echo <-
#'   get.adjust(stata = haven::read_dta('data/ARICmaster_121820.dta'),
#'              mods  = list('data/visit-five/echovars.csv'))
#' }
get.adjust <- function(stata, mods, filts = NULL) {

  args    <- purrr::map_df(mods,  readr::read_csv)
  stack   <- paste0('`', args[[1]], '` = ', args[[2]], collapse = ',')
  mutates <- paste0('dplyr::mutate(stata,', stack, ')')
  stata   <- parse(text = mutates) %>% eval()

  if (!is.null(filts)) {
    stata <- get.filtered(stata, filts)
  }

  stata <- dplyr::select(stata, dplyr::all_of(args[[1]]))
  labelled::var_label(stata) <- paste0(args[[3]])
  return(stata)
}


#' Filter Data given series of filter commands
#'
#' This function is a helper function to take a series of filter instructions
#' (from the \code{filters.csv} file) and uses them to filter a tibble of study
#' data.
#'
#' This command writes an intermediate csv file to the working directory, named
#' \code{filres.csv}. This intermediate file stores the data about how many
#' subjects were filtered by each condition. This data is later used to create
#' consort diagrams.
#'
#' @param data The Tibble of Data from Study to be Filtered
#' @param filts The list of file path(s) to \code{filters.csv}. Most likely a
#'   list of length 1.
#'
#' @return returns a filtered tibble of study data. Does not do any error checking.
#' @export
#' @importFrom magrittr %>%
#' @examples
#' \dontrun{
#'
#' filts <- '~/proteomics/filters.csv'
#' data  <- haven::read_dta('~/proteomics/studydata.dta')
#'
#' filtered.data <- get.filtered(data, filts)
#'
#' OR
#'
#' data <- haven::read_dta('~/proteomics/studydata.dta') %>%
#'     get.filtered('~/proteomics/filters.csv')
#'
#' }
#'
get.filtered <- function(data, filts) {
  filt  <- purrr::map_df(filts, readr::read_csv)
  calls <- paste0('`', filt[[1]], '` ', filt[[2]], ' ', filt[[3]])
  filters <- paste0('dplyr::filter(data, ', calls, ')')

  sizes <- purrr::map(filters, ~ dim(data)[[1]] - dim(parse(text = .x) %>% eval())[[1]])
  excls <- paste0(calls, ' -- ', sizes)

  calls   <- paste0(calls, collapse = ',')
  filters <- paste0('dplyr::filter(data, ', calls, ')')
  stata   <- parse(text = filters) %>% eval()
  excls   <- c(paste0("Initial -- ", dim(data)[[1]]),
               excls, paste0("Final -- ", dim(stata)[[1]])) %>% as.data.frame()

  readr::write_csv(excls, 'filres.csv', col_names = FALSE)

  return(stata)
}

#' Data Joining, Imputation, and Scaling
#'
#' @param master Master study data tibble from [haven::read_dta()]
#' @param soma.data Parsed SomaLogic proteomic data tibble from [get.soma()]
#' @param adjust List of Adjustment variables. Must match \code{adjusted.csv}
#'   files, but should not include id variable.
#'
#' @return Tibble of joined and cleaned proteomics and study data.
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' \dontrun{
#' labels <- get.labels('data/proteinMapping.csv', sid,
#'                      uniprot.full.name,
#'                      entrezgenesymbol,
#'                      flag2 == 0)
#'
#' soma.data.fifth <- get.soma('data/somaV5.txt', labels = labels)
#'
#' fifth.visit.study <-
#'   get.adjust(stata = haven::read_dta('data/ARICmaster_121820.dta'),
#'              mods  = list('data/visit-five/outcomes.csv',
#'                           'data/visit-five/adjusted.csv'),
#'              filts = list('data/visit-five/filters.csv'))
#'
#' adj <- readr::read_csv('data/visit-five/adjusted.csv')[[1]][-1]
#'
#' visit.data <- get.data(fifth.visit.study, soma.data.fifth, adj)
#' }
get.data <- function(master, soma.data, adjust) {

  imputed <- master %>%
    dplyr::select(dplyr::all_of(adjust)) %>%
    mice::mice(printFlag = FALSE)

  missing <- imputed$nmis[imputed$nmis != 0]
  howmany <-
    tibble::tibble(term = names(missing),
                   miss = missing) %>%
    dplyr::transmute(val = paste0(term, ' -- ', miss))
  readr::write_csv(howmany, 'impres.csv', col_names = F)

  imputed <- imputed %>%
    mice::complete() %>%
    tibble::as_tibble()

  outcomes <- master %>%
    dplyr::select(!dplyr::all_of(adjust)) %>%
    dplyr::bind_cols(imputed)

  data <- dplyr::inner_join(outcomes, soma.data, by = "id")

  return(data)
}


#' Scale Data
#'
#' This function scales data after tables are made. Scales numeric adjustment
#' and proteomic variables. Done separately to maintain scale and values for
#' table 1 and figure generation. Numeric adjustment variables must be
#' specified. See [numeric.id()] for more details.
#'
#' @param data Tibble of Data to be scaled
#' @param adjust List of NUMERIC ADJUSTMENT VARS ONLY -- output of
#'   [numeric.id()].
#' @param proteins List of Proteins to Scale
#'
#' @return Tibble where all numeric adjustment variables and protein variables
#'   are scaled to mean 0 and sd 1.
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' \dontrun{
#' labels <- get.labels('data/proteinMapping.csv', sid,
#'                      uniprot.full.name,
#'                      entrezgenesymbol,
#'                      flag2 == 0)
#'
#' soma.data.fifth <- get.soma('data/somaV5.txt', labels)
#'
#' fifth.visit.study <-
#'   get.adjust(stata = haven::read_dta('data/ARICmaster_121820.dta'),
#'              mods  = list('data/visit-five/outcomes.csv',
#'                           'data/visit-five/adjusted.csv'),
#'              filts = list('data/visit-five/filters.csv'))
#'
#' adj <- readr::read_csv('data/visit-five/adjusted.csv')[[1]][-1]
#'
#' visit.data <- get.data(fifth.visit.study, soma.data.fifth, adj)
#'
#' num.vars <- numeric.id(visit.data, adj)
#' visit.data.scaled <- get.scaled(visit.data, num.vars, labels$term)
#' }
get.scaled <- function(data, adjust, proteins = NULL) {

  data %>%
    dplyr::mutate_at(dplyr::all_of(c(adjust, proteins)), ~ scale(.x)[,1]) %>%
    return()
}


#' Get Master Visit Data File with Proteomics
#'
#' This is the outward facing wrapper function to process and load data for
#' analysis. This function will read SomaLogic and Study Data, Impute missing
#' adjusting covariates, filter data for exclusions
#'
#' @param soma.data SomaLogic Proteomic Raw Data from [get.soma()]
#' @param path.mods Char vector path to folder structure containing new
#'   adjustment variables, filtering conditions, and outcome variables. This
#'   folder must contain csv files formatted as detailed in [get.adjust()]. This
#'   function expects that adjustment variables will be in \code{adjusted.csv},
#'   outcome variables in \code{outcomes.csv} and filters in \code{filters.csv}.
#' @param master Tibble of master data file with study data. Output of
#'   [haven::read_dta()] function.
#'
#' @return Tibble of study data merged with Proteomics Data, labelled,
#'   and filtered.
#'
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' \dontrun{
#' labels <- get.labels('data/proteinMapping.csv', sid,
#'                      uniprot.full.name,
#'                      entrezgenesymbol,
#'                      flag2 == 0)
#' aric.master <- haven::read_dta('data/ARICmaster_121820.dta')
#' soma.data <- get.soma('data/somaV5.txt', labels)
#'
#' fifth.visit <- get.visit(soma.data = soma.data,
#'                          path.mods = 'data/visit-five',
#'                          labels = labels,
#'                          master = aric.master))
#' }

get.visit <- function(soma.data, path.mods, master) {

  visit.data <- master %>% dplyr::filter(id %in% soma.data$id)
  adtl.data <- get.adjust(stata = visit.data,
                          mods  = list(file.path(path.mods, 'outcomes.csv'),
                                       file.path(path.mods, 'adjusted.csv')),
                          filts = list(file.path(path.mods, 'filters.csv')))

  adjust <- readr::read_csv(file.path(path.mods, 'adjusted.csv'))[[1]]
  visit.data <- get.data(adtl.data, soma.data, adjust)
  return(visit.data)
}


#' Get Numeric Variables (for scaling)
#'
#' Numeric Values will likely need to be scaled for survival analysis. This
#' function identifies which variables are numeric, separating them from
#' categorical factor variables. used in [get.visit()] and [add.data()]
#'
#' @param data Tibble of data containing variables to be considered
#' @param vars Vector of all variables to be considered
#'
#' @return Vector containing column names of numeric variables
#'
#' @importFrom magrittr %>%
#'
#' @export
#'
#' @examples
#' \dontrun{
#' new.vars <- readr::read_csv('data/visit-five/echovars.csv')[[1]]
#' the.data <- haven::read_dta('data/ARICmaster_121820.dta')
#' num.vars <- numeric.id(vars = new.vars, data = the.data)
#' }
numeric.id <- function(data, vars) {
  type <- name <- NULL
  num <- data %>%
    dplyr::select(dplyr::all_of(vars)) %>%
    dplyr::summarise_all(class) %>%
    tidyr::gather('name', 'type') %>%
    dplyr::filter(type == 'numeric') %>%
    dplyr::pull(name)

  return(num)
}

#' Augment Data Variables
#'
#' This function serves to add variables to the master tibble dataframe.
#' Variables are joined in on the id column.
#'
#' @details See [get.adjust()] for csv formatting details.
#'
#' @param new.path Path to Stata .dta file containing new variables
#' @param old.data Tibble of old study data to merge with
#' @param scaled Boolean indicating whether new data should be scaled
#' @param mods List of paths to formatted csv files creating new adjustment
#'   variables. Format of csv is detailed below. This parameter is
#'   required. See [get.adjust()] for formatting.
#' @param filts List of path (singular) to formatted csv of the data filters to
#'   be applied to the dataset. See [get.adjust()] for formatting.
#'
#' @return Tibble with new variables joined in.
#'
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' \dontrun{
#' fifth.visit.adtl <- add.data(new.path = 'data/ARIC_proteomics_addlvar_01062021.dta',
#'                              old.data = haven::read_dta('data/ARICmaster_121820.dta'),
#'                              mods = list('data/visit-five/adtlvars.csv'))
#' }
add.data <- function(new.path, old.data, scaled = F, mods, filts = NULL) {

  new.vars <- readr::read_csv(mods)[[1]]
  new.data <- haven::read_dta(new.path) %>%
    get.adjust(mods = mods, filts = filts) %>%
    dplyr::right_join(old.data, by = "id")

  new.data <-
    dplyr::relocate(new.data, new.vars[-1],
                    .before = grep("SeqId", names(new.data), value = T)[1])

  if (scaled) {
    num.var <- numeric.id(data = new.data, vars = new.vars)
    new.data <- get.scaled(new.data, num.var)
  }

  return(new.data)
}

