library(ggplot2)

consort <- function(log.path, mod.path, soma.data, visit.data, labels, fupt, hfdiag) {

  small.height <- 1.75
  large.height <- 7
  small.width <- 3.5
  large.width <- 2
  spacing <- 0.75

  num.prot <- glue::glue('Protein quantification by modified aptamer assay
                         {prot} Proteins measured in a total of {samp} samples',
                         prot = dim(soma.data)[2] - 1,
                         samp = dim(soma.data)[1])

  num.qc <- glue::glue('{prot} Proteins, {samp} Samples Passed Quality Control',
                       prot = dim(labels)[1],
                       samp = dim(soma.data)[1])

  num.excl <-
    glue::glue('Exclusions:\n',
               paste0('  ',
                      readr::read_csv(file.path(log.path, 'filres.csv'),
                                      col_names = FALSE) %>%
                  dplyr::pull(),
                  collapse = '\n'),
               '\nImputed:\n',
               paste0('  ',
                      readr::read_csv(file.path(log.path, 'impres.csv'),
                                      col_names = FALSE) %>%
                      dplyr::pull(),
                      collapse = '\n'),
               .trim = F)
  num.data <-
    glue::glue('Visit Dataset:\nN = {samp}, n = {evnt}, Median Follow Up = {medn} years',
               samp = dim(visit.data)[1],
               evnt = sum(visit.data[[hfdiag]]),
               medn = round(median(visit.data[[fupt]]) / 365.25, 2))

  text <- c(num.prot, num.qc, num.excl, num.data)

  boxes <- tibble::tibble(b = c("s", "l", "s", "s")) %>%
    dplyr::mutate(x1 = dplyr::if_else(b == "s", 1, 0),
                  x2 = dplyr::if_else(b == "s", x1 + small.width, x1 + large.width),
                  y2 = dplyr::if_else(b == "s", small.height, large.height),
                  hj = dplyr::if_else(b == "s",0.5, 0),
                  sp = spacing,
                  sp = cumsum(sp),
                  y3 = cumsum(y2) + sp,
                  y1 = c(0, y3[-length(b)]),
                  y2 = y1 + y2,
                  lx = dplyr::if_else(b == "s",x1+(x2-x1)/2, x1+(x2-x1)/6),
                  ly = y1+(y2-y1)/2,
                  tx = text[length(b):1]) %>%
    dplyr::select(x1, x2, y1, y2, tx, hj, lx, ly)

  arrow <- dplyr::select(boxes, y1, y2) %>%
    dplyr::slice(-which(boxes$x1 == 0)) %>%
    dplyr::mutate(y2 = c(0, y2[-length(y2)]),
                  x1 = small.width/2 + 1) %>%
    dplyr::slice(-1)

  p <-
    ggplot() +
    geom_rect(data=boxes,
              mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2),
              fill=rgb(255,255,247, maxColorValue = 255),
              color="black",
              alpha=1) +
    geom_text(data=boxes,
              aes(x=lx, y=ly, label=tx, hjust=hj),
              size=4) +
    geom_segment(data = arrow,
                 aes(x = x1, y = y1 - 0.1, xend = x1, yend = y2 + 0.1),
                 arrow = arrow(length=unit(0.25, "cm"))) +
    geom_segment(aes(x = small.width/2 + 1, xend = large.width + 0.05,
                     y    = boxes[which(boxes$x1 == 0), "ly"][[1]],
                     yend = boxes[which(boxes$x1 == 0), "ly"][[1]]),
                 arrow = arrow(length = unit(0.25, "cm"))) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())

    return(p)
}




